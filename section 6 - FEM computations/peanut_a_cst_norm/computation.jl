using Pkg: Pkg
Pkg.activate(".")
Pkg.instantiate()

using DelimitedFiles, SparseArrays, LinearAlgebra
using LinearMaps, IterativeSolvers

function load_matrix(filename::String)
    data = readdlm(filename, ' ', Float64, '\n')
    rows, cols = Int.(data[:, 1]), Int.(data[:, 2])
    vals = Complex{Float64}.(data[:, 3], data[:, 4])

    nb_dof = maximum(rows)
    return sparse(rows, cols, vals, nb_dof, nb_dof)
end

# function norm_inv(A::AbstractMatrix)
#     N = size(A, 1)
#     vals, _, _, _ = svdsolve(A, 1, :SR; maxiter = N)
#     return 1 / minimum(vals)
# end

function norm_inv(A::AbstractMatrix)
    F = factorize(A)
    F_ad = adjoint(F)

    function H_inv_yx(y, x)
        w = zeros(eltype(x), size(x))
        ldiv!(w, F, x)
        return ldiv!(y, F_ad, w)
    end

    H_inv_map = LinearMap{eltype(A)}(H_inv_yx, size(A, 1); ismutating = true,
        ishermitian = true)

    return √powm(H_inv_map)[1]
end

function load_resonances(filename::String)
    data = readdlm(filename, ' ', Float64, '\n')
    return Complex{Float64}.(data[:, 1], data[:, 2]), Complex{Float64}.(data[:, 3], data[:, 4])
end

function make_sample(k1::Vector{<:Number}, k2::Vector{<:Number}, nb_inter::Integer)
    kr = (real.(k1) + real.(k2)) / 2
    km = append!([0.0], (kr[1:(end-1)] .+ kr[2:end]) / 2)
    km = append!(km, km[end] + (km[end] - km[end-1]))

    sample_k = []

    for i in 1:(length(km)-1)
        v = range(km[i], km[i+1]; length = nb_inter + 2)[2:(end-1)]
        append!(sample_k, v)
        append!(sample_k, km[i+1])
    end

    h = minimum(sample_k[2:end] - sample_k[1:(end-1)])
    for i in 1:15
        v = kr[i] .+ [-h / 4, 0, h / 4]
        append!(sample_k, v)
        append!(sample_k, real.([k1[i], k2[i]]))
    end

    return sort!(sample_k)
end

function main()
    eMin = "-1.2"
    eMax = "-1.1"
    for dGL in 1:8
        filename = "disk_eVar_$(eMin)_$(eMax)_$dGL"
        nb_inter = 4 # even

        open(filename, "w") do io
            return writedlm(io, [])
        end

        k_res_1, k_res_2 = load_resonances("disk_eVar_$(eMin)_$(eMax)_pla")
        sample_k = make_sample(k_res_1, k_res_2, nb_inter)
        display(sample_k)

        run(`make`)

        for k in sample_k
            run(`./exec-x86_64-linux-g++-10-Release $eMin $eMax $k $dGL`)
            A = load_matrix("tmp_matrix")
            ni = norm_inv(A)
            open(filename, "a") do io
                return writedlm(io, [k real(ni) imag(ni)])
            end
        end
    end

    return rm("tmp_matrix")
end

@time main()
