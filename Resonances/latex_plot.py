""" Function for pretty save """


def set_size(width=370.38374, frac_width=1, frac_height=None):
    # Width of figure (in pts)
    fig_width = width * frac_width
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    inv_golden_ratio = (5**0.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width * inches_per_pt
    # Figure height in inches
    if frac_height is None:
        fig_height_in = fig_width_in * inv_golden_ratio
    else:
        fig_height_in = fig_width_in * frac_height

    return (fig_width_in, fig_height_in)


def set_rcParams(font_size=10, line_width=1):
    parameters = {
        # Use LaTeX to write all text
        "text.usetex": True,
        "font.family": "serif",
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": font_size,
        "font.size": font_size,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": font_size,
        "xtick.major.width": 0.5,
        "xtick.labelsize": font_size - 2,
        "ytick.major.width": 0.5,
        "ytick.labelsize": font_size - 2,
        # Line properties
        "lines.linewidth": line_width,
        "axes.linewidth": 0.5,
        "grid.linewidth": 0.25,
    }
    return parameters
