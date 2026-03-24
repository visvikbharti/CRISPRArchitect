## CRISPRArchitect v0.1.0 — Initial Release

A computational framework for predicting HDR gene conversion outcomes and optimizing multi-site genome editing strategies.

### Modules
- **ConversionSim** — Monte Carlo HDR gene conversion tract simulator
- **MOSAIC** — Multi-locus editing strategy optimizer (HDR/base/prime editing)
- **cssDNA-TopoPred** — ssDNA secondary structure analyzer
- **ChromBridge** — 3D chromatin distance & translocation risk predictor
- **LoopSim** — Cohesin loop extrusion simulator

### Features
- Streamlit web interface with dark/light mode
- Ensembl REST API integration (real gene fetching, auto donor design)
- CLI interface for command-line usage
- Docker container for reproducible deployment
- 30/30 unit tests passing

### Validation
- ConversionSim validated against 4 published datasets (Iyer 2022, Chauhan 2023, Elliott 1998, Paquet 2016)
- MOSAIC benchmarked against 14 published genome editing studies (71.4% concordance)

### Quick Start
```bash
git clone https://github.com/visvikbharti/CRISPRArchitect.git
cd CRISPRArchitect
pip install numpy scipy matplotlib
streamlit run webapp/app.py
```

### Citation
Bharti V and Chakraborty D. CRISPRArchitect (2026). https://github.com/visvikbharti/CRISPRArchitect
