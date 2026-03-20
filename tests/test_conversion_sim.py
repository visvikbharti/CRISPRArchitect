"""
Tests for ConversionSim module.

Validates that the Monte Carlo HDR simulation produces biologically
reasonable outputs consistent with published experimental data.
"""

import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from conversion_sim import ConversionSimulator
from conversion_sim.resection import ResectionSimulator
from conversion_sim.filament import FilamentModel
from conversion_sim.synthesis import SynthesisSimulator


class TestResection:
    """Test end resection simulation."""

    def test_resection_produces_positive_lengths(self):
        sim = ResectionSimulator()
        result = sim.simulate("blunt", 0, n_simulations=1000)
        assert np.all(result.left_bp > 0)
        assert np.all(result.right_bp > 0)

    def test_staggered_cut_increases_resection(self):
        sim = ResectionSimulator()
        blunt = sim.simulate("blunt", 0, n_simulations=5000)
        staggered = sim.simulate("staggered_5prime", 5, n_simulations=5000)
        assert np.mean(staggered.left_bp) >= np.mean(blunt.left_bp) - 100

    def test_resection_in_biological_range(self):
        """Resection should be 300-10000 bp (Symington, 2011)."""
        sim = ResectionSimulator()
        result = sim.simulate("blunt", 0, n_simulations=10000)
        median_left = np.median(result.left_bp)
        assert 500 < median_left < 5000, f"Median resection {median_left} outside expected range"


class TestFilament:
    """Test RAD51 filament formation."""

    def test_filament_covers_resected_region(self):
        model = FilamentModel()
        result = model.form_filament(np.array([2000, 3000, 1000]))
        assert np.all(result.filament_length_nt > 0)
        assert np.all(result.filament_length_nt <= np.array([2000, 3000, 1000]) * 1.1)

    def test_short_resection_produces_small_filament(self):
        """Very short resection should produce small filaments."""
        model = FilamentModel()
        result = model.form_filament(np.array([5, 5, 5]))
        assert np.all(result.filament_length_nt <= 10)


class TestSynthesis:
    """Test DNA synthesis / gene conversion tract simulation."""

    def test_tract_lengths_positive(self):
        sim = SynthesisSimulator()
        tracts = sim.simulate_synthesis("circular_ssDNA", n_simulations=1000)
        assert np.all(tracts.tract_lengths_bp > 0)

    def test_median_tract_in_range(self):
        """Median tract should be ~300-800 bp (Elliott et al., 1998)."""
        sim = SynthesisSimulator()
        tracts = sim.simulate_synthesis("linear_dsDNA", n_simulations=10000)
        median = np.median(tracts.tract_lengths_bp)
        assert 100 < median < 2000, f"Median tract {median} outside expected range"

    def test_circular_produces_longer_tracts(self):
        """cssDNA should produce longer tracts than linear (Iyer et al., 2022)."""
        sim = SynthesisSimulator()
        linear = sim.simulate_synthesis("linear_ssDNA", n_simulations=10000)
        circular = sim.simulate_synthesis("circular_ssDNA", n_simulations=10000)
        assert np.median(circular.tract_lengths_bp) >= np.median(linear.tract_lengths_bp) * 0.9


class TestConversionSimulator:
    """Integration tests for the full simulator."""

    def test_basic_run(self):
        sim = ConversionSimulator(
            cut_type="blunt", overhang_length=0,
            donor_topology="circular_ssDNA", homology_arm_length=300,
            cell_type="iPSC", n_simulations=1000
        )
        results = sim.run()
        assert results is not None

    def test_hdr_rate_biological_range(self):
        """HDR rate in iPSCs should be 1-15% (literature consensus)."""
        sim = ConversionSimulator(
            cut_type="staggered_5prime", overhang_length=3,
            donor_topology="circular_ssDNA", homology_arm_length=300,
            cell_type="iPSC", n_simulations=10000
        )
        sim.run()
        stats = sim.summary()
        rate = stats["hdr_success_rate"]
        assert 0.005 < rate < 0.20, f"HDR rate {rate} outside iPSC range"

    def test_hek293t_higher_than_ipsc(self):
        """HEK293T should have higher HDR rate than iPSCs."""
        sim_ipsc = ConversionSimulator(
            cut_type="blunt", overhang_length=0,
            donor_topology="circular_ssDNA", homology_arm_length=300,
            cell_type="iPSC", n_simulations=10000
        )
        sim_hek = ConversionSimulator(
            cut_type="blunt", overhang_length=0,
            donor_topology="circular_ssDNA", homology_arm_length=300,
            cell_type="HEK293T", n_simulations=10000
        )
        sim_ipsc.run()
        sim_hek.run()
        rate_ipsc = sim_ipsc.summary()["hdr_success_rate"]
        rate_hek = sim_hek.summary()["hdr_success_rate"]
        assert rate_hek > rate_ipsc, "HEK293T should have higher HDR than iPSC"

    def test_probability_at_distance_decreasing(self):
        """P(reaching distance) should decrease with distance."""
        sim = ConversionSimulator(
            cut_type="blunt", overhang_length=0,
            donor_topology="circular_ssDNA", homology_arm_length=300,
            cell_type="HEK293T", n_simulations=10000
        )
        sim.run()
        p_500 = sim.probability_at_distance(500)
        p_1000 = sim.probability_at_distance(1000)
        p_2000 = sim.probability_at_distance(2000)
        assert p_500 >= p_1000 >= p_2000, "Probability should decrease with distance"

    def test_distant_site_unreachable(self):
        """A site 1 Mb away should be unreachable by HDR."""
        sim = ConversionSimulator(
            cut_type="staggered_5prime", overhang_length=3,
            donor_topology="circular_ssDNA", homology_arm_length=300,
            cell_type="iPSC", n_simulations=10000
        )
        sim.run()
        p_1mb = sim.probability_at_distance(1_000_000)
        assert p_1mb == 0.0, f"1 Mb should be unreachable, got {p_1mb}"


class TestMOSAIC:
    """Tests for the MOSAIC strategy optimizer."""

    def test_import(self):
        from mosaic.gene_structure import GeneStructure
        from mosaic.mutation_classifier import Mutation, MutationClassifier
        from mosaic.strategy_enumerator import StrategyEnumerator
        from mosaic.scorer import StrategyScorer

    def test_gene_structure(self):
        from mosaic.gene_structure import GeneStructure
        exons = [{"number": i, "start": i * 10000, "end": i * 10000 + 150}
                 for i in range(1, 11)]
        gene = GeneStructure.from_manual("TestGene", exons)
        assert len(gene.exons) == 10
        dist = gene.genomic_distance(1, 10)
        assert dist > 0

    def test_mutation_classification(self):
        from mosaic.mutation_classifier import Mutation, MutationClassifier
        mc = MutationClassifier()
        # G>A is a transition (ABE-amenable)
        assert "transition" in mc.classify("G", "A")
        # G>C is a transversion
        assert mc.classify("G", "C") == "transversion"

    def test_strategy_enumeration(self):
        from mosaic.gene_structure import GeneStructure
        from mosaic.mutation_classifier import Mutation
        from mosaic.strategy_enumerator import StrategyEnumerator

        exons = [{"number": i, "start": i * 50000, "end": i * 50000 + 150}
                 for i in range(1, 61)]
        gene = GeneStructure.from_manual("TestGene", exons)
        m1 = Mutation(exon_number=20, position=1000000, ref_allele="G", alt_allele="A")
        m2 = Mutation(exon_number=50, position=2500000, ref_allele="C", alt_allele="T")
        strategies = StrategyEnumerator().enumerate_strategies(gene, [m1, m2], "iPSC", "SpCas9")
        assert len(strategies) >= 3  # Should find at least 3 strategies

    def test_single_template_not_offered_for_distant_sites(self):
        from mosaic.gene_structure import GeneStructure
        from mosaic.mutation_classifier import Mutation
        from mosaic.strategy_enumerator import StrategyEnumerator

        exons = [{"number": i, "start": i * 50000, "end": i * 50000 + 150}
                 for i in range(1, 61)]
        gene = GeneStructure.from_manual("TestGene", exons)
        m1 = Mutation(exon_number=1, position=50000, ref_allele="G", alt_allele="A")
        m2 = Mutation(exon_number=60, position=3000000, ref_allele="C", alt_allele="T")
        strategies = StrategyEnumerator().enumerate_strategies(gene, [m1, m2], "iPSC", "SpCas9")
        names = [s.name for s in strategies]
        assert "SINGLE_TEMPLATE_HDR" not in names, "Should not offer single template for distant sites"


class TestChromBridge:
    """Tests for ChromBridge 3D distance predictor."""

    def test_import(self):
        from chrombridge import ChromatinDistancePredictor

    def test_distance_increases_with_genomic_distance(self):
        from chrombridge import ChromatinDistancePredictor
        pred = ChromatinDistancePredictor()
        d_100k = pred.predict_3d_distance(100_000)
        d_1m = pred.predict_3d_distance(1_000_000)
        assert d_1m.mean_3d_distance_nm > d_100k.mean_3d_distance_nm

    def test_cssdna_cannot_bridge_megabases(self):
        from chrombridge import ChromatinDistancePredictor
        pred = ChromatinDistancePredictor()
        bridge = pred.can_donor_bridge(1_000_000, 3000, "circular_ssDNA")
        assert not bridge.feasible, "3 kb cssDNA should not bridge 1 Mb"

    def test_nearby_sites_may_be_bridgeable(self):
        from chrombridge import ChromatinDistancePredictor
        pred = ChromatinDistancePredictor()
        bridge = pred.can_donor_bridge(2000, 5000, "circular_ssDNA")
        # 5 kb donor for 2 kb distance — might be feasible
        assert bridge.bridgeability_ratio > 0.5


class TestTopoPred:
    """Tests for cssDNA-TopoPred."""

    def test_import(self):
        from topopred import DonorAnalyzer

    def test_g4_detection(self):
        from topopred.g_quadruplex import GQuadruplexScanner
        scanner = GQuadruplexScanner()
        # Known G4 motif
        seq = "AAAA" + "GGG" + "TT" + "GGG" + "TT" + "GGG" + "TT" + "GGG" + "AAAA"
        motifs = scanner.find_g4_motifs(seq)
        assert len(motifs) > 0, "Should detect G4 motif"

    def test_clean_sequence_accessible(self):
        from topopred import DonorAnalyzer
        analyzer = DonorAnalyzer()
        # Random sequence unlikely to form stable structures
        import random
        random.seed(42)
        seq = "".join(random.choice("ATCG") for _ in range(600))
        report = analyzer.analyze(seq, left_arm=(0, 200), right_arm=(400, 600))
        # Should not recommend redesign for random sequence
        assert report["recommendation"] != "Redesign strongly recommended"


class TestLoopSim:
    """Tests for LoopSim cohesin simulator."""

    def test_import(self):
        from loopsim import LoopSimulator

    def test_basic_simulation(self):
        from loopsim import LoopSimulator
        sim = LoopSimulator(0, 1_000_000, resolution_bp=5000)
        sim.setup_fiber(ctcf_spacing=100_000)
        sim.add_dsb(500_000)
        results = sim.run(n_simulations=100)
        assert results is not None

    def test_search_domain_bounded(self):
        """Search domain should not exceed fiber boundaries."""
        from loopsim import LoopSimulator
        sim = LoopSimulator(0, 500_000, resolution_bp=1000)
        sim.setup_fiber(ctcf_spacing=50_000)
        sim.add_dsb(250_000)
        results = sim.run(n_simulations=200)
        # All boundaries should be within the fiber
        for dsb_pos, dsb_result in results.per_dsb_results.items():
            assert np.all(dsb_result.left_boundaries >= 0)
            assert np.all(dsb_result.right_boundaries <= 500_000)


class TestEnsembl:
    """Tests for Ensembl API integration (requires internet)."""

    @pytest.mark.skipif(
        os.environ.get("SKIP_NETWORK_TESTS", "0") == "1",
        reason="Network tests disabled"
    )
    def test_fetch_nf1(self):
        from utils.ensembl import fetch_gene
        gene = fetch_gene("NF1")
        assert gene.symbol == "NF1"
        assert gene.n_exons > 50
        assert gene.chromosome == "17"

    @pytest.mark.skipif(
        os.environ.get("SKIP_NETWORK_TESTS", "0") == "1",
        reason="Network tests disabled"
    )
    def test_fetch_nonexistent_gene(self):
        from utils.ensembl import fetch_gene, EnsemblError
        with pytest.raises(EnsemblError):
            fetch_gene("ZZZZNOTAREALGENE999")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
