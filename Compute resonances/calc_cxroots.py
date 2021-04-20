""" Compute resonances using the cxroots library (contour integration techniques)

    Authors: Zoïs Moitier, Camille Carvalho
            Karlsruhe Institute of Technology, Germany
            University of California, Merced

    Last modified: 20/04/2021
"""

from sys import argv

import matplotlib.pyplot as plt
import numpy as np
from cxroots import AnnulusSector, Circle
from scipy.special import h1vp, hankel1, iv, ivp

## Entries ##
ε = float(argv[1]) #For example -1.1 + 1e-2 * 1j 
η = np.sqrt(-ε)
print(f"η = {η}")
c = η + 1 / η


## Internal functions ##
def rootsAnnSec(m, rMin, rMax, aMin, aMax):
    f0 = lambda k: ivp(m, η * k) * hankel1(m, k) / η + iv(m, η * k) * h1vp(m, k)
    f1 = (
        lambda k: ivp(m, η * k, 2) * hankel1(m, k)
        + c * ivp(m, η * k) * h1vp(m, k)
        + iv(m, η * k) * h1vp(m, k, 2)
    )

    A = AnnulusSector(center=0.0, radii=(rMin, rMax), phiRange=(aMin, aMax))
    z = A.roots(f0, df=f1)
    return z.roots


def writeFile(myFile, m, z):
    if np.size(z, 0):
        for i in range(np.size(z, 0)):
            myFile.write(f"{m} {z[i].real} {z[i].imag}\n")


def calcInt():
    plaTrue = ε > -1.0

    if plaTrue:
        Int = open(f"eps_{ε}_int", "w")
        Pla = open(f"eps_{ε}_pla", "w")
    else:
        Int = open(f"eps_{ε}_int", "w")

    for m in range(65):
        print(f"m = {m}")

        f0 = lambda k: ivp(m, η * k) * hankel1(m, k) / η + iv(m, η * k) * h1vp(m, k)
        f1 = (
            lambda k: ivp(m, η * k, 2) * hankel1(m, k)
            + c * ivp(m, η * k) * h1vp(m, k)
            + iv(m, η * k) * h1vp(m, k, 2)
        )

        t = np.linspace(0.2, 65.0, num=1024)
        k = 1j * t
        rf = np.real(f0(k))

        ind = np.where(rf[1:] * rf[:-1] < 0.0)[0]
        roots = np.zeros(np.shape(ind), dtype=complex)
        for a, i in enumerate(ind):
            C = Circle(center=1j * (t[i] + t[i + 1]) / 2.0, radius=(t[i + 1] - t[i]))
            z = C.roots(f0, df=f1)
            roots[a] = z.roots[0]

        if plaTrue:
            if m:
                writeFile(Int, m, roots[1:])
                writeFile(Pla, m, roots[[0]])
            else:
                writeFile(Int, m, roots)
        else:
            writeFile(Int, m, roots)

    if plaTrue:
        Int.close()
        Pla.close()
    else:
        Int.close()

calcInt()

def calcResPla():
    if ε < -1.0:
        Pla = open(f"eps_{ε}_pla", "w")

        angle = -np.pi / 4.0
        for m in range(1, 65):
            r = max(0.1, 0.9 * np.sqrt(1.0 - η ** (-2)) * m - 1.0)
            R = max(2.0, 1.1 * np.sqrt(1.0 - η ** (-2)) * m + 1.0)
            a = min(angle, -1e-3)

            z = rootsAnnSec(m, r, R, a, 1e-3)
            writeFile(Pla, m, z)
            angle = np.angle(z[0])

        Pla.close()

calcResPla()


def calcResOut():
    Out = open(f"eps_{ε}_out", "w")

    rMin = 0.2
    rMax = 5.0
    aMin = -np.pi + 0.01
    aMax = 0.0

    for m in range(33, 65):
        print(f"m = {m}")

        z = rootsAnnSec(m, rMin, rMax, aMin, aMax)
        writeFile(Out, m, z)

        if m > 3:
            zMod = np.abs(z)
            zArg = np.angle(z)
            rMin = max(0.2, np.amin(zMod) * 0.75)
            rMax = max(rMax, np.amax(zMod) + 3.0)
            aMin = min(aMin, (-np.pi + np.amin(zArg)) / 2.0)
            aMax = np.amax(zArg) / 2.0

    Out.close()


calcResOut()


def calc_cx_pla():
    with open(f"eps_{ε}_pla", "w") as file:
        rMin, rMax = 0.1, 0.5
        aMin = -np.pi / 4
        for m in range(1, 65):
            z = rootsAnnSec(m, rMin, rMax, aMin, 1e-3)[0]
            file.write(f"{m} {z.real} {z.imag}\n")

            rMin = abs(z)
            rMax = abs(z) * (m + 1) / m + 1
            aMin = min(2.5 * np.angle(z), -1e-3)
            print(m, rMin, rMax, aMin)


calc_cx_pla()


def rewriteSave():
    Int = np.loadtxt(f"eps_{ε}_int")
    Pla = np.loadtxt(f"eps_{ε}_pla")
    Out = np.loadtxt(f"eps_{ε}_out")

    ind = np.argsort(Out[:, 1])[::-1]
    out2 = Out[ind]
    rep = out2[:, 1] > -1e-3

    np.savez(f"eps_{ε}.npz", inner=Int, plasmon=Pla, outer=out2[rep])

rewriteSave()


def rewriteSave_pla():
    Pla = np.loadtxt(f"eps_{ε}_pla")
    np.savez(f"eps_{ε}.npz", plasmon=Pla)


#rewriteSave_pla()
