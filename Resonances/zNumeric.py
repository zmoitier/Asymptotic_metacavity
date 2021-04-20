import numba
import numpy as np
from scipy import integrate as si
from scipy import special as sp


@numba.vectorize([numba.float64(numba.complex128), numba.float32(numba.complex64)])
def abs2(x):
    return x.real ** 2 + x.imag ** 2


def ε_to_η(ε):
    return np.sqrt(-ε)


def det_M0(η, m):
    return lambda z: (
        sp.ivp(m, η * z) * sp.hankel1(m, z) / η + sp.iv(m, η * z) * sp.h1vp(m, z)
    )


def det_M1(η, m):
    c = η + 1 / η
    return lambda z: (
        sp.ivp(m, η * z, 2) * sp.hankel1(m, z)
        + c * sp.ivp(m, η * z) * sp.h1vp(m, z)
        + sp.iv(m, η * z) * sp.h1vp(m, z, 2)
    )


def asy_res_m(η, m):
    return m * m * (1 - η ** (-2))


def res_field_m(η, m, k):
    iv, h1 = sp.iv(m, η * k), sp.hankel1(m, k)
    return (lambda r: sp.iv(m, η * k * r) / iv, lambda r: sp.hankel1(m, k * r) / h1)


def res_field_m_far(η, m, k):
    iv, h1 = sp.iv(m, η * k), sp.hankel1(m, k)
    return (
        lambda r: sp.iv(m, η * k * r) / iv,
        lambda r: np.exp(1j * k * r) * sp.hankel1e(m, k * r) / h1,
    )


def asy_field_m(η, m):
    return (
        lambda r: np.exp(η * m * (r - 1)),
        lambda r: np.exp(-m * (r - 1) / η),
    )


def asyH1(m, k, r):
    z = k * r
    ω = z - m * np.pi / 2 - np.pi / 4
    return np.sqrt(2 / (np.pi * z)) * np.exp(1j * ω) / sp.hankel1(m, k)


def coeff_w_m(η, m, k):
    M = np.array(
        [[sp.iv(m, η * k), -sp.hankel1(m, k)], [sp.ivp(m, η * k) / η, sp.h1vp(m, k)]]
    )
    V = np.array([[sp.jv(m, k)], [-sp.jvp(m, k)]])
    S = np.linalg.solve(M, V)
    return (S[0, 0], S[1, 0])


def total_field_m(η, m, k):
    α, β = coeff_w_m(η, m, k)
    return (
        lambda r: α * sp.iv(m, η * k * r),
        lambda r: β * sp.hankel1(m, k * r) + sp.jv(m, k * r),
    )


def scattered_field_m(η, m, k):
    α, β = coeff_w_m(η, m, k)
    return (
        lambda r: α * sp.iv(m, η * k * r) - sp.jv(m, k * r),
        lambda r: β * sp.hankel1(m, k * r),
    )


def field_r(field, R):
    w = np.zeros(R.shape, dtype=complex)
    Int, Out = np.where(np.less_equal(R, 1)), np.where(np.greater(R, 1))
    w[Int] += field[0](R[Int])
    w[Out] += field[1](R[Out])
    return w


def res_field_xy(η, m, k, X, Y):
    R, T = np.hypot(X, Y), np.arctan2(Y, X)
    return field_r(res_field_m(η, m, k), R) * np.exp(1j * m * T)


def normalization(u, R):
    Int = R <= 1
    ind = np.unravel_index(np.argmax(np.abs(u[Int])), u[Int].shape)
    return u / u[Int][ind]


def field_xy(η, k, field, X, Y):
    R, T = np.hypot(X, Y), np.arctan2(Y, X)

    M = 32  # 1 + np.where(sp.jv(np.arange(100), np.amax(R) * k) > 1e-8)[0][-1]
    print(f"M = {M}")

    U = np.zeros(R.shape, dtype=complex)
    for m in range(M, 0, -1):
        if m % 2:
            U += complex(0, 2) * field_r(field(η, m, k), R) * np.sin(m * T)
        else:
            U += 2 * field_r(field(η, m, k), R) * np.cos(m * T)
    U += field_r(field(η, 0, k), R)
    return U


def normL2_radial2(field, T):
    f = lambda r: abs2(field[0](r)) * r
    g = lambda r: abs2(field[1](r)) * r
    return si.quad(f, 0, 1)[0] + si.quad(g, 1, T)[0]


def N_ερ(η, k, ρ):
    M = 32  # 1 + np.where(sp.jv(np.arange(100), ρ * k) > 1e-8)[0][-1]
    result = 0
    for m in range(M, 0, -1):
        result += 2 * normL2_radial2(scattered_field_m(η, m, k), ρ)
    result += normL2_radial2(scattered_field_m(η, 0, k), ρ)
    return np.sqrt(2 * result) / ρ


def res_load_mq(ε, m, q):
    data = np.load(f"../data/eps_{ε}.npz")

    Int = data["inner"]
    Pla = data["plasmon"]
    Out = data["outer"]

    vm = np.concatenate((Int[:, 0], Pla[:, 0], Out[:, 0]), axis=0).astype(int)
    re = np.concatenate((Int[:, 1], Pla[:, 1], Out[:, 1]), axis=0)
    im = np.concatenate((Int[:, 2], Pla[:, 2], Out[:, 2]), axis=0)

    ind = np.where(vm == abs(m))[0]

    return complex(re[ind[q]], im[ind[q]])


def sample_geo(kMin, kMax, Nlin, val, Ngeo, δ):
    k = np.linspace(kMin, kMax, num=Nlin)

    sample = list(k)
    vδ = np.geomspace(δ, 1, num=Ngeo, endpoint=False)
    I = np.searchsorted(k, val)
    for i, v in enumerate(val[val <= kMax]):
        dk = min(v - k[I[i] - 2], k[I[i] + 1] - v)

        sample.extend(list(v - dk * vδ))
        sample.append(v)
        sample.extend(list(v + dk * vδ))

    print(f"N = {len(sample)}")
    return np.sort(sample)


def calc_response(ε, vk, ρ):
    η = ε_to_η(ε)
    vn = np.zeros(vk.shape)
    for i, k in enumerate(vk):
        vn[i] = N_ερ(η, k, ρ)
    return vn
