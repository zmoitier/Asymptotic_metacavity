"""Plot mode from FEM"""
import sys

import meshio
import numpy as np
from matplotlib import pyplot as plt


def ri_to_c(filename: str) -> np.array:
    val = np.loadtxt(filename)
    return val[:, 2] + 1j * val[:, 3]


def print_val() -> None:
    val0 = ri_to_c("data/eigs_s0")
    val1 = ri_to_c("data/eigs_s1")
    print(np.abs(val0 - val1).max())

    return None


def plot_mode(filename: str) -> None:
    reader = meshio.read(f"{filename}")

    x = reader.points

    quad = reader.cells[0].data
    tri = np.zeros((2 * np.size(quad, 0), 3), dtype=int)
    for q in range(np.size(quad, 0)):
        tri[2 * q, :] = [quad[q, 0], quad[q, 1], quad[q, 2]]
        tri[2 * q + 1, :] = [quad[q, 0], quad[q, 2], quad[q, 3]]

    term = list(reader.point_data.keys())
    print("-" * len(f"Real part .... ← {term[0]}"))
    print(f"Real part .... ← {term[0]}")
    print(f"Imaginary part ← {term[1]}")
    print(f"Modulus ...... ← {term[2]}")

    ur = reader.point_data[term[0]]
    ui = reader.point_data[term[1]]
    um = reader.point_data[term[2]]

    plt.subplot(2, 2, 1)
    urM = np.amax(np.abs(ur))
    plt.tricontourf(x[:, 0], x[:, 1], tri, ur, 128, cmap="RdBu_r")
    plt.axis("equal")
    plt.clim(-urM, urM)
    plt.colorbar()
    plt.title("Real part", fontsize=15)

    plt.subplot(2, 2, 2)
    uiM = np.amax(np.abs(ui))
    plt.tricontourf(x[:, 0], x[:, 1], tri, ui, 128, cmap="RdBu_r")
    plt.axis("equal")
    plt.clim(-uiM, uiM)
    plt.colorbar()
    plt.title("Imaginary part", fontsize=15)

    plt.subplot(2, 2, 3)
    umM = np.amax(um)
    plt.tricontourf(x[:, 0], x[:, 1], tri, um, 128)
    plt.axis("equal")
    plt.clim(0.0, umM)
    plt.colorbar()
    plt.title("Modulus", fontsize=15)

    plt.subplot(2, 2, 4)
    ua = np.arctan2(ui, ur)
    plt.tricontourf(x[:, 0], x[:, 1], tri, ua, 128, cmap="twilight_shifted")
    plt.axis("equal")
    cbar = plt.colorbar(ticks=[-np.pi, -np.pi / 2.0, 0, np.pi / 2.0, np.pi])
    cbar.ax.set_yticklabels(
        [r"$-\pi$", r"$-\frac{\pi}{2}$", r"$0$", r"$\frac{\pi}{2}$", r"$\pi$"]
    )
    plt.title("Argument", fontsize=15)

    plt.show()

    return None


if __name__ == "__main__":
    # print_val()

    filename = sys.argv[1]
    plot_mode(filename)
