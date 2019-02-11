module HelmholtzEigenBase

using ArnoldiFormat,
ArnoldiMethod,
CoordinateSystems,
BoundaryConditions,
DifferentialOperators,
LinearAlgebra,
SparseArrays

export eig_kl,
eig_cf


"""
    k, ψ =  eig_kl(N, dx, ρ, k, bcs; nk=1, ka=0, kb=0, coordinate_system=:cart, h=1)

number of lattice sites `N`, lattice spacing `dx`, density function `ρ`, boundary conditions `bcs` (see `DifferentialOperators.ezbc`),
number of frequencies `nk`, Floquet/Bloch wavenumbers `ka`,`kb`.
"""
function eig_kl(ρ, N, ∂Ω, dx, coordinate_system::Type{T}, bc, bl, k;
                nk::Int=1, ka::Number=0, kb::Number=0, lattice_index::Int=1,
                disp_opt::Bool = true) where {Th, T<:CoordinateSystem}

    (∇², fs, s), defs = laplacian(N,coordinate_system,∂Ω,bc,bl, k, ka, kb, lattice_index)
    A = sum(∇²)
    A = (A + SparseMatrixCSC(transpose(A)))/2
    B = s*spdiagm(0=>-ρ[:])
    M = shift_and_invert(A, B, k^2, diag_inv_B=true, issymmetric=true)
    decomp, history = partialschur(M, nev=nk)

    disp_opt ? println(history) : nothing
    @assert history.converged history

    # this is left to user to do on own because partialeigen is currently type unstable
    if !issymmetric(M)
        @warn "differential operator not symmetric, pass decomp (third output argument) to ArnoldiMethod.partialeigen"
    end

    k² = k^2 .+ 1 ./decomp.eigenvalues
    ψ = Array(decomp.Q)

    for i ∈ 1:nk
        𝒩² = sum((ψ[:,i].^2).*diag(B))*(isinf(dx[1]) ? 1 : dx[1])*(isinf(dx[2]) ? 1 : dx[2])
        @views ψ[:,i] = ψ[:,i]/sqrt(𝒩²)/exp(complex(0,angle(ψ[end÷2-1,i])))
    end
    return sqrt.(k²), ψ, decomp
end


"""
    η, u =  eig_cf(N, dx, ρ, k, bcs; nk=1, ka=0, kb=0, coordinate_system=:cart, h=1)

`F` and `ρ` must be of the same size.
"""
function eig_cf(ρ, F, N, ∂Ω, dx, coordinate_system::Type{T}, bc, bl, k;
                η::Number=0, ncf::Int=1,
                ka::Number=0, kb::Number=0, lattice_index::Int=1,
                disp_opt::Bool = true
                ) where {Th, T<:CoordinateSystem}

    (∇², fs, s), defs = laplacian(N,Cartesian,∂Ω,bc,bl, k, ka, kb, lattice_index)
    A = sum(∇²)
    A = sum(∇²) + s*spdiagm(0 => ρ[:]*k^2)
    A = (A + SparseMatrixCSC(transpose(A)))/2
    B = s*spdiagm(0 => -F[:]*k^2)
    M = shift_and_invert(A, B, η, diag_inv_B=true, issymmetric=true)
    decomp, history = partialschur(M, nev=ncf)

    disp_opt ? println(history) : nothing
    @assert history.converged history

    # this is left to user to do on own because partialeigen is currently type unstable
    if !issymmetric(M)
        @warn "differential operator not symmetric, pass decomp (third output argument) to ArnoldiMethod.partialeigen"
    end

    η = η .+ 1 ./decomp.eigenvalues
    u = Array(decomp.Q)

    for i ∈ 1:ncf
        𝒩² = sum((u[:,i].^2).*diag(B))*(isinf(dx[1]) ? 1 : dx[1])*(isinf(dx[2]) ? 1 : dx[2])
        @views u[:,i] = u[:,i]/sqrt(𝒩²)/exp(complex(0,angle(u[end÷2-1,i])))
    end
    return η, u, decomp
end


end # module
