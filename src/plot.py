""" Utils for ploting

    Author: Zoïs Moitier
            Karlsruhe Institute of Technology, Germany
"""
from math import sqrt


def set_size(width=6.4, frac_width=1, frac_height=None):
    # Width of figure (in pts)
    fig_width = width * frac_width

    # Golden ratio to set aesthetic figure height
    inv_golden_ratio = (sqrt(5) - 1) / 2

    # Figure height in inches
    if frac_height is None:
        fig_height = fig_width * inv_golden_ratio
    else:
        fig_height = fig_width * frac_height

    return (fig_width, fig_height)


def set_rcParams(font_size=10):
    parameters = {
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": font_size,
        "font.size": font_size,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": font_size,
        "xtick.major.width": 0.5,
        "xtick.labelsize": font_size - 2,
        "ytick.major.width": 0.5,
        "ytick.labelsize": font_size - 2,
    }
    return parameters
