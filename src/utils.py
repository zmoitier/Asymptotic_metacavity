""" Utils for ploting

    Author: Zoïs Moitier
            Karlsruhe Institute of Technology, Germany
"""


def is_number(s):
    """ Returns True is string is a number. """
    try:
        float(s)
        return True
    except ValueError:
        return False
