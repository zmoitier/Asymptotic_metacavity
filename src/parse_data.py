""" Utils for ploting

    Author: Zoïs Moitier
            Karlsruhe Institute of Technology, Germany
"""
import numpy as np


def parse_data_k(data, type, M):
    values = data[type]
    m = values[:, 0].astype(int)
    ind_m = np.where(np.equal(m, M))[0]
    return (m, values[ind_m, 1] + 1j * values[ind_m, 2])


def parse_data_λ(data, type, M):
    values = data[type]
    m = values[:, 0].astype(int)
    ind_m = np.where(np.less_equal(m, M))[0]
    return (m, np.power(values[ind_m, 1] + 1j * values[ind_m, 2], 2))
