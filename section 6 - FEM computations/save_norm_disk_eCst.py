"""Save graph reponse."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullLocator

import src


def main():
    ε = -1.1

    plt.rcParams.update(src.set_rcParams(font_size=8, line_width=1))
    fig, ax = plt.subplots(
        figsize=src.set_size(frac_width=0.49), constrained_layout=True
    )

    data_fem = np.loadtxt("../data/disk_eCst_-1.1_8")
    k_fem = data_fem[:, 0]
    r_fem = data_fem[:, 1]

    data = np.load(f"../data/eps_{ε}.npz")
    vk = src.sample_geo(0.01, 5, 512, data["plasmon"][:, 1], 3, 0.25)
    vn = [src.norm_inverse(ε, k, 1, 32) for k in vk]

    for k_res in data["plasmon"][:, 1]:
        ax.semilogy([k_res, k_res], [1, 1e10], color="#b0b0b0", lw=0.25, alpha=1)
    ax.semilogy(vk, vn, "C0", label="analytic")
    ax.semilogy(k_fem, r_fem, "C1--", label="FEM")

    ax.set_xlim(0, k_fem[-1])
    ax.set_xlabel("$k$")

    ax.set_ylim(1, 1e10)
    ax.set_ylabel("revolvent norm")

    # ax.grid(visible=True, which="major", axis="y")
    ax.legend(loc=2)

    # plt.show()

    ext = "pdf"
    fig.savefig(f"../images/disk_eCst_norm.{ext}", transparent=True, format=ext)
    plt.close(fig=fig)

    return None


if __name__ == "__main__":
    main()
