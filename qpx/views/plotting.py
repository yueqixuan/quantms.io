"""Plot helpers for QPX views — lazy matplotlib import."""

from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import matplotlib.figure


def _get_plt():
    """Lazy matplotlib import with non-interactive backend."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def bar_chart(
    df,
    x: str,
    y: str,
    title: str = "",
    xlabel: str = "",
    ylabel: str = "",
    figsize: tuple = (10, 6),
    rotation: int = 45,
) -> "matplotlib.figure.Figure":
    """Simple bar chart from a DataFrame."""
    plt = _get_plt()
    fig, ax = plt.subplots(figsize=figsize)
    ax.bar(df[x].astype(str), df[y])
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=rotation)
    fig.tight_layout()
    return fig


def grouped_bar_chart(
    df,
    x: str,
    y_cols: list[str],
    title: str = "",
    xlabel: str = "",
    ylabel: str = "",
    figsize: tuple = (12, 6),
) -> "matplotlib.figure.Figure":
    """Grouped bar chart with multiple y-columns."""
    import numpy as np

    plt = _get_plt()
    fig, ax = plt.subplots(figsize=figsize)
    x_vals = df[x].astype(str)
    n_groups = len(y_cols)
    width = 0.8 / n_groups
    for i, col in enumerate(y_cols):
        positions = np.arange(len(x_vals)) + i * width
        ax.bar(positions, df[col], width, label=col)
    ax.set_xticks(np.arange(len(x_vals)) + width * (n_groups - 1) / 2)
    ax.set_xticklabels(x_vals, rotation=45, ha="right")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend()
    fig.tight_layout()
    return fig
