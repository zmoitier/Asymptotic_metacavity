Λ_term[0] = (η - 1) * (η + 1) / η**2

ϑ_term[0] = 1

P_term[0] = 1

Q_term[0] = 1

Λ_term[1] = -((η - 1) ** 2) * (η + 1) ** 2 * κ(s) / η**3

ϑ_term[1] = (Λ[1] * η**3 + η**4 * κ(s) - 2 * η**2 * κ(s) + κ(s)) / (
    2 * η * (η - 1) * (η + 1)
)

P_term[1] = (
    σ * (Λ[1] * η**4 - η**2 * σ * κ(s) + σ * κ(s)) / (2 * η * (η - 1) * (η + 1))
)

Q_term[1] = -η * σ * (Λ[1] - η**2 * σ * κ(s) + σ * κ(s)) / (2 * (η - 1) * (η + 1))

Λ_term[2] = -(
    -Λ[1] ** 2 * η**8
    + 2 * η**12 * κ(s) ** 2
    + 2 * I * η**11 * Derivative(κ(s), s)
    - 5 * η**10 * κ(s) ** 2
    - 4 * I * η**9 * Derivative(κ(s), s)
    + 6 * η**8 * κ(s) ** 2
    + 2 * I * η**7 * Derivative(κ(s), s)
    - 6 * η**6 * κ(s) ** 2
    - 2 * I * η**5 * Derivative(κ(s), s)
    + 6 * η**4 * κ(s) ** 2
    + 4 * I * η**3 * Derivative(κ(s), s)
    - 5 * η**2 * κ(s) ** 2
    - 2 * I * η * Derivative(κ(s), s)
    + 2 * κ(s) ** 2
) / (4 * η**6 * (η - 1) * (η + 1))

ϑ_term[2] = (
    -Λ[1] ** 2 * η**8
    + 4 * Λ[2] * η**8
    - 4 * Λ[2] * η**6
    + 2 * η**12 * κ(s) ** 2
    + 2 * I * η**11 * Derivative(κ(s), s)
    - 5 * η**10 * κ(s) ** 2
    - 4 * I * η**9 * Derivative(κ(s), s)
    + 6 * η**8 * κ(s) ** 2
    + 2 * I * η**7 * Derivative(κ(s), s)
    - 6 * η**6 * κ(s) ** 2
    - 2 * I * η**5 * Derivative(κ(s), s)
    + 6 * η**4 * κ(s) ** 2
    + 4 * I * η**3 * Derivative(κ(s), s)
    - 5 * η**2 * κ(s) ** 2
    - 2 * I * η * Derivative(κ(s), s)
    + 2 * κ(s) ** 2
) / (8 * η**4 * (η - 1) ** 2 * (η + 1) ** 2)

P_term[2] = (
    σ
    * (
        3 * Λ[1] ** 2 * η**10 * σ
        - 3 * Λ[1] ** 2 * η**9
        - 6 * Λ[1] * η**8 * σ**2 * κ(s)
        - 6 * Λ[1] * η**7 * σ * κ(s)
        + 6 * Λ[1] * η**6 * σ**2 * κ(s)
        + 6 * Λ[1] * η**5 * σ * κ(s)
        + 12 * Λ[2] * η**9
        - 12 * Λ[2] * η**7
        + 6 * η**11 * κ(s) ** 2
        + 6 * I * η**10 * Derivative(κ(s), s)
        - 12 * η**9 * κ(s) ** 2
        - 6 * η**8 * σ * κ(s) ** 2
        - 18 * I * η**8 * Derivative(κ(s), s)
        + 12 * η**7 * σ**2 * κ(s) ** 2
        + 6 * I * η**7 * σ * Derivative(κ(s), s)
        + 12 * η**7 * κ(s) ** 2
        + 3 * η**6 * σ**3 * κ(s) ** 2
        + 4 * I * η**6 * σ**2 * Derivative(κ(s), s)
        + 12 * η**6 * σ * κ(s) ** 2
        + 18 * I * η**6 * Derivative(κ(s), s)
        - 28 * η**5 * σ**2 * κ(s) ** 2
        - 18 * I * η**5 * σ * Derivative(κ(s), s)
        - 12 * η**5 * κ(s) ** 2
        - 6 * η**4 * σ**3 * κ(s) ** 2
        - 8 * I * η**4 * σ**2 * Derivative(κ(s), s)
        - 6 * I * η**4 * Derivative(κ(s), s)
        + 20 * η**3 * σ**2 * κ(s) ** 2
        + 18 * I * η**3 * σ * Derivative(κ(s), s)
        + 6 * η**3 * κ(s) ** 2
        + 3 * η**2 * σ**3 * κ(s) ** 2
        + 4 * I * η**2 * σ**2 * Derivative(κ(s), s)
        - 12 * η**2 * σ * κ(s) ** 2
        - 4 * η * σ**2 * κ(s) ** 2
        - 6 * I * η * σ * Derivative(κ(s), s)
        + 6 * σ * κ(s) ** 2
    )
    / (24 * η**4 * (η - 1) ** 2 * (η + 1) ** 2)
)

Q_term[2] = (
    σ
    * (
        3 * Λ[1] ** 2 * η**6
        + 3 * Λ[1] ** 2 * η**5 * σ
        + 6 * Λ[1] * η**8 * σ * κ(s)
        - 6 * Λ[1] * η**7 * σ**2 * κ(s)
        - 6 * Λ[1] * η**6 * σ * κ(s)
        + 6 * Λ[1] * η**5 * σ**2 * κ(s)
        - 12 * Λ[2] * η**6
        + 12 * Λ[2] * η**4
        + 6 * η**11 * σ * κ(s) ** 2
        + 4 * η**10 * σ**2 * κ(s) ** 2
        + 6 * I * η**10 * σ * Derivative(κ(s), s)
        + 3 * η**9 * σ**3 * κ(s) ** 2
        + 4 * I * η**9 * σ**2 * Derivative(κ(s), s)
        - 12 * η**9 * σ * κ(s) ** 2
        - 20 * η**8 * σ**2 * κ(s) ** 2
        - 18 * I * η**8 * σ * Derivative(κ(s), s)
        - 6 * η**8 * κ(s) ** 2
        - 6 * η**7 * σ**3 * κ(s) ** 2
        - 8 * I * η**7 * σ**2 * Derivative(κ(s), s)
        - 6 * I * η**7 * Derivative(κ(s), s)
        + 28 * η**6 * σ**2 * κ(s) ** 2
        + 18 * I * η**6 * σ * Derivative(κ(s), s)
        + 12 * η**6 * κ(s) ** 2
        + 3 * η**5 * σ**3 * κ(s) ** 2
        + 4 * I * η**5 * σ**2 * Derivative(κ(s), s)
        + 12 * η**5 * σ * κ(s) ** 2
        + 18 * I * η**5 * Derivative(κ(s), s)
        - 12 * η**4 * σ**2 * κ(s) ** 2
        - 6 * I * η**4 * σ * Derivative(κ(s), s)
        - 12 * η**4 * κ(s) ** 2
        - 6 * η**3 * σ * κ(s) ** 2
        - 18 * I * η**3 * Derivative(κ(s), s)
        + 12 * η**2 * κ(s) ** 2
        + 6 * I * η * Derivative(κ(s), s)
        - 6 * κ(s) ** 2
    )
    / (24 * η**3 * (η - 1) ** 2 * (η + 1) ** 2)
)
