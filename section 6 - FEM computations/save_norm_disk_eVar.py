"""Save graph reponse."""

import matplotlib.pyplot as plt
import numpy as np

import src


def main():
    εMin, εMax = -1.2, -1.1

    plt.rcParams.update(src.set_rcParams(font_size=8, line_width=1))
    fig, ax = plt.subplots(
        figsize=src.set_size(frac_width=0.49), constrained_layout=True
    )

    data_pla = np.loadtxt("../data/disk_eVar_-1.2_-1.1_pla")

    data_fem = np.loadtxt("../data/disk_eVar_-1.2_-1.1_8")
    k_fem = data_fem[:, 0]
    r_fem = data_fem[:, 1]

    for k_res in data_pla[:, 2]:
        ax.semilogy([k_res, k_res], [1, 1e10], color="#b0b0b0", lw=0.25, alpha=1)
    ax.semilogy(k_fem, r_fem, "C1--")

    ax.set_xlim(0, k_fem[-1])
    ax.set_xlabel("$k$")

    ax.set_ylim(1e3, 1e10)
    ax.set_ylabel("revolvent norm")

    # ax.grid(visible=True, which="major", axis="y")

    # plt.show()

    ext = "pdf"
    fig.savefig(f"../images/disk_eVar_norm.{ext}", transparent=True, format=ext)
    plt.close(fig=fig)

    return None


if __name__ == "__main__":
    main()
