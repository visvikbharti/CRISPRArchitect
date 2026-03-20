"""
CRISPRArchitect Webapp — Custom CSS and Styling Helpers
========================================================

Fully compatible with BOTH Streamlit light and dark themes.
Users toggle via Settings (hamburger menu) > Theme > Light/Dark.

All custom elements use CSS custom properties (variables) that
automatically adapt to whichever theme Streamlit is using.
"""

import streamlit as st


# ---------------------------------------------------------------------------
# Colour constants (for Python-side use in charts etc.)
# ---------------------------------------------------------------------------

PRIMARY = "#2196F3"
SUCCESS = "#4CAF50"
WARNING = "#FF9800"
DANGER = "#F44336"
INFO = "#00BCD4"
SURFACE = "#F5F7FA"
TEXT_DARK = "#1A1A2E"
TEXT_MUTED = "#6C757D"

GRADE_COLOURS = {
    "A": "#4CAF50",
    "B": "#8BC34A",
    "C": "#FF9800",
    "D": "#FF5722",
    "F": "#F44336",
}

SCORE_COLOUR_MAP = {
    "efficiency": "#2196F3",
    "safety": "#4CAF50",
    "time": "#FF9800",
    "cost": "#9C27B0",
    "overall": "#1A1A2E",
}


# ---------------------------------------------------------------------------
# Inject global CSS — theme-aware (light + dark)
# ---------------------------------------------------------------------------

def inject_global_css() -> None:
    """Inject theme-aware CSS that works in both light and dark mode.

    Streamlit uses CSS custom properties like:
      --background-color, --secondary-background-color, --text-color
    We hook into those so our custom elements adapt automatically.

    The user can toggle themes via the hamburger menu > Settings > Theme.
    """
    st.markdown(
        """
        <style>
        /* ================================================
           BASE STYLES (adapt to Streamlit's theme vars)
           ================================================ */

        .block-container {
            max-width: 1100px;
            padding-top: 2rem;
        }

        /* ---- Card containers ---- */
        .ca-card {
            background: var(--secondary-background-color, #FFFFFF);
            border: 1px solid rgba(128, 128, 128, 0.2);
            border-radius: 10px;
            padding: 1.2rem 1.4rem;
            margin-bottom: 1rem;
            box-shadow: 0 1px 4px rgba(0, 0, 0, 0.08);
        }
        .ca-card h3 {
            margin-top: 0;
            color: var(--text-color, #1A1A2E);
            font-size: 1.1rem;
        }

        /* ---- Grade badge ---- */
        .grade-badge {
            display: inline-flex;
            align-items: center;
            justify-content: center;
            width: 48px;
            height: 48px;
            border-radius: 50%;
            color: white;
            font-weight: 700;
            font-size: 1.4rem;
            box-shadow: 0 2px 6px rgba(0, 0, 0, 0.2);
        }

        /* ---- Section headers ---- */
        .section-header {
            font-size: 1.25rem;
            font-weight: 600;
            color: #2196F3;
            border-bottom: 2px solid #2196F3;
            padding-bottom: 0.35rem;
            margin-top: 1.5rem;
            margin-bottom: 0.8rem;
        }

        /* ---- Small help labels ---- */
        .field-help {
            font-size: 0.82rem;
            color: var(--text-color, #6C757D);
            opacity: 0.7;
            margin-bottom: 0.2rem;
        }

        /* ---- Strategy row highlights ---- */
        .strategy-best {
            background: rgba(76, 175, 80, 0.15);
            border-left: 4px solid #4CAF50;
            padding: 0.5rem 0.8rem;
            border-radius: 6px;
            margin-bottom: 0.5rem;
            color: var(--text-color);
        }
        .strategy-alt {
            background: rgba(255, 152, 0, 0.12);
            border-left: 4px solid #FF9800;
            padding: 0.5rem 0.8rem;
            border-radius: 6px;
            margin-bottom: 0.5rem;
            color: var(--text-color);
        }
        .strategy-risky {
            background: rgba(244, 67, 54, 0.12);
            border-left: 4px solid #F44336;
            padding: 0.5rem 0.8rem;
            border-radius: 6px;
            margin-bottom: 0.5rem;
            color: var(--text-color);
        }

        /* ---- Metric card row ---- */
        .metric-row {
            display: flex;
            gap: 1rem;
            flex-wrap: wrap;
            margin-bottom: 1rem;
        }
        .metric-item {
            background: var(--secondary-background-color, #F5F7FA);
            border: 1px solid rgba(128, 128, 128, 0.15);
            border-radius: 8px;
            padding: 0.8rem 1rem;
            flex: 1;
            min-width: 140px;
            text-align: center;
        }
        .metric-item .label {
            font-size: 0.8rem;
            color: var(--text-color, #6C757D);
            opacity: 0.65;
            text-transform: uppercase;
            letter-spacing: 0.05em;
        }
        .metric-item .value {
            font-size: 1.5rem;
            font-weight: 700;
            color: var(--text-color, #1A1A2E);
        }

        /* ================================================
           SIDEBAR — let Streamlit handle colors natively.
           Only customize font sizes and spacing.
           ================================================ */

        [data-testid="stSidebar"] h1 {
            font-size: 1.3rem;
            color: #2196F3 !important;
        }

        </style>
        """,
        unsafe_allow_html=True,
    )


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def section_header(title: str) -> None:
    """Render a styled section header."""
    st.markdown(
        f'<div class="section-header">{title}</div>',
        unsafe_allow_html=True,
    )


def card_start() -> None:
    """Open a card container (must pair with ``card_end``)."""
    st.markdown('<div class="ca-card">', unsafe_allow_html=True)


def card_end() -> None:
    """Close a card container."""
    st.markdown("</div>", unsafe_allow_html=True)


def grade_badge(grade: str) -> str:
    """Return an HTML snippet for a coloured letter-grade badge."""
    colour = GRADE_COLOURS.get(grade.upper(), "#757575")
    return (
        f'<span class="grade-badge" style="background:{colour};">'
        f"{grade.upper()}</span>"
    )


def render_grade(grade: str, label: str = "Safety Grade") -> None:
    """Display a grade badge with an optional label using ``st.markdown``."""
    st.markdown(
        f"{grade_badge(grade)}&nbsp;&nbsp;<strong>{label}</strong>",
        unsafe_allow_html=True,
    )


def score_colour(value: float) -> str:
    """Return a CSS colour string for a 0-1 score."""
    if value >= 0.70:
        return SUCCESS
    elif value >= 0.40:
        return WARNING
    else:
        return DANGER


def format_bp(bp) -> str:
    """Format a base-pair count with appropriate units."""
    bp = float(bp)
    if bp >= 1_000_000:
        return f"{bp / 1_000_000:.2f} Mb"
    elif bp >= 1_000:
        return f"{bp / 1_000:.1f} kb"
    else:
        return f"{int(bp):,} bp"


def format_pct(value: float, decimals: int = 1) -> str:
    """Format a 0-1 fraction as a percentage string."""
    return f"{value * 100:.{decimals}f}%"


def metric_row_html(items) -> str:
    """Build an HTML metric row from (label, value) tuples."""
    inner = ""
    for label, value in items:
        inner += (
            f'<div class="metric-item">'
            f'<div class="label">{label}</div>'
            f'<div class="value">{value}</div>'
            f"</div>"
        )
    return f'<div class="metric-row">{inner}</div>'


def strategy_css_class(rank: int, safety_score: float) -> str:
    """Return the CSS class name for a strategy card."""
    if rank == 1:
        return "strategy-best"
    elif safety_score < 0.35:
        return "strategy-risky"
    else:
        return "strategy-alt"
