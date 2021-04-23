import sympy as sy


def apply_var(fct_iter, *var):
    return list(map(lambda f: f(*var), fct_iter))


def metric_expan(z, h, N):
    if N == 1:
        return (1, 1)
    else:
        g = 1 + z * h
        ig1 = sy.series(1 / g, x=h, n=N)
        return (g, ig1.removeO())


def inv_ε_expan(η_iter, z, h, N):
    ηiη0 = 0
    for q, η in enumerate(η_iter):
        ηiη0 += η * (z * h) ** q / (η_iter[0] * sy.factorial(q))
    εiε0 = sy.series(ηiη0 ** 2, x=h, n=N)
    ε0iε = sy.series(1 / εiε0, x=h, n=N)
    return (-ε0iε / η_iter[0] ** 2).expand().removeO()


def formal_expan(fn, h, N):
    expan = fn[0]
    for n in range(1, N):
        expan += fn[n] * h ** n
    return expan


def op_disk(iε, dz, z, h, N):
    f = sy.symbols("f", cls=sy.Function)
    r, ir = metric_expan(z, h, N)
    tmp = (iε * r * dz(f(z)) + sy.O(h ** N)).expand()
    L = (-ir * dz(tmp) + iε * ir * ir * f(z) + sy.O(h ** N)).expand()
    return (f, L.removeO())


def op_omega(iε, κ, ds, dz, s, z, h, N):
    f, t = sy.symbols("f t", cls=sy.Function)
    g, ig = metric_expan(κ(s) * z, h, N)
    tmp_s = (iε * ig * ds(f(s, z), t(s)) + sy.O(h ** N)).expand()
    tmp_z = (iε * g * dz(f(s, z)) + sy.O(h ** N)).expand()
    L = (-ig * (ds(tmp_s, t(s)) + dz(tmp_z)) + sy.O(h ** N)).expand()
    return (f, t, L.removeO())


def list_solve_exp(d, a, z):
    """Compute solution to -P''-2aP' = zʲ for j in {0,d}"""
    M = sy.zeros(d + 1, d + 1)
    M[0, 0] = -2 * a
    for j in range(1, d + 1):
        M[j, j] = -2 * a
        M[j - 1, j] = -j
    N = M.inv()
    list_exp = []
    for j in range(d, -1, -1):
        list_exp.append(sum([N[i, j] * z ** (i + 1) / (i + 1) for i in range(j + 1)]))
    return list_exp


def solve_exp(a, z, S):
    S_coeffs = sy.Poly(S, z).all_coeffs()
    d = len(S_coeffs) - 1
    list_exp = list_solve_exp(d, a, z)
    return (d + 1, sum([S_coeffs[i] * list_exp[i] for i in range(d + 1)]).expand())


def factor_term(expr):
    f_expr = 0
    for t in expr.expand().args:
        f_expr += t.factor()
    return f_expr


def real_imag_parts(expr):
    re = expr.subs(sy.I, 0).expand()
    im = (expr - re).subs(sy.I, 1).expand()
    return (re, im)


def save_Λ(file_name, Λ, h, N):
    with open(f"save_{file_name}.py", "w") as file:
        for n in range(N):
            term = Λ.coeff(h, n)
            term_real = term.subs(sy.I, 0).factor()
            term_imag = (term - term_real).factor().subs(sy.I, 1)
            file.write(f"term_real[{n}] = {term_real}\n\n")
            file.write(f"term_imag[{n}] = {term_imag}\n\n")


def save_omega(file_name, Λ, ϑ, P, Q, h, N):
    with open(f"save_{file_name}.py", "w") as file:
        for n in range(N):
            file.write(f"Λ_term[{n}] = {Λ.coeff(h, n).factor()}\n\n")
            file.write(f"ϑ_term[{n}] = {ϑ.coeff(h, n).factor()}\n\n")
            file.write(f"P_term[{n}] = {P.coeff(h, n).factor()}\n\n")
            file.write(f"Q_term[{n}] = {Q.coeff(h, n).factor()}\n\n")
