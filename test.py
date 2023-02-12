import matplotlib.pyplot as plt
import numpy as np
from scipy.special import ai_zeros

from src import parse_data_k

ε = -1.1
η = np.sqrt(-ε)
data = np.load(f"data/eps_{ε}.npz")

nb = 59
k = np.zeros(nb, dtype=complex)
for m in range(1, nb + 1):
    k_inn = parse_data_k(data, "inner", m)
    k[m - 1] = k_inn[0]

aj = -ai_zeros(1)[0]

m = np.arange(1, nb + 1)
h = 2 / m
asy1 = 1j * m / η
asy2 = 1j * m * (1 + aj * h ** (2 / 3) / 2) / η
asy3 = 1j * m * (1 + aj * h ** (2 / 3) / 2 - h ** 2 / (2 * η * np.sqrt(η ** 2 + 1))) / η


def slope(x, y):
    return (y[1:] - y[:-1]) / (x[1:] - x[:-1])


def rate(x, y):
    return slope(np.log(x), np.log(y))


plt.plot(m[1:], rate(m, np.abs(asy1 / k - 1)))
plt.plot(m[1:], rate(m, np.abs(asy2 / k - 1)))
plt.plot(m[1:], rate(m, np.abs(asy3 / k - 1)))

plt.show()
