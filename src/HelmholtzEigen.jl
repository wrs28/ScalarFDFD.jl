module HelmholtzEigen

using ArnoldiHelper,
Arpack,
ArnoldiMethod,
DifferentialOperators,
LinearAlgebra,
SparseArrays

export eig_kl,
eig_cf


"""
    k, ψ =  eig_kl(N, dx, ρ, k, bcs; nk=1, ka=0, kb=0, coordinate_system=:cart, h=1)
    k, ψ =  eig_kl(dx, ρ, k, bcs; nk=1, ka=0, kb=0, coordinate_system=:cart, h=1)

number of lattice sites `N`, lattice spacing `dx`, density function `ρ`, boundary conditions `bcs` (see `DifferentialOperators.ezbc`),
number of frequencies `nk`, Floquet/Bloch wavenumbers `ka`,`kb`.
if `N` if omitted, it is inferred from density `ρ`
"""
function eig_kl(N, dx, ρ, k::Number, bcs::Tuple{T,U};
                nk::Int=1, ka::Number=0, kb::Number=0,
                coordinate_system::Symbol=:cart, h=ones(N...)) where U<:BoundaryCondition where T<:BoundaryCondition

    ∇², S = laplacian(N, dx, bcs; ka=ka, kb=kb, coordinate_system=coordinate_system, h=h)
    decomp, history = partialschur(shift_and_invert(∇², spdiagm(0=>-ρ[:]), k^2, diag_inv_B=true), nev=nk)
    history.converged ? nothing : @warn "incomplete convergence: only $(h.nconverged) of $(h.nev) evecs converged"
    k², ψ = k^2 .+ 1 ./decomp.eigenvalues, decomp.Q
    for i ∈ 1:nk
        𝒩² = sum(abs2.(ψ[:,i].*ρ[:]))*prod(dx)
        @views ψ[:,i] = ψ[:,i]/sqrt(𝒩²)
    end
    return sqrt.(k²), ψ
end
function eig_kl(dx, ρ, k::Number, bcs::Tuple{T,U};
                nk::Int=1, ka::Number=0, kb::Number=0,
                coordinate_system::Symbol=:cart, h=1) where U<:BoundaryCondition where T<:BoundaryCondition
    return eig_kl([size(ρ)...], dx, ρ, k, bcs; nk=nk, ka=ka, kb=kb, coordinate_system=coordinate_system, h=h)
end


"""
    η, u =  eig_cf(N, dx, ρ, k, bcs; nk=1, ka=0, kb=0, coordinate_system=:cart, h=1)
    η, u =  eig_cf(dx, ρ, k, bcs; nk=1, ka=0, kb=0, coordinate_system=:cart, h=1)

if `N` if omitted, it is inferred from `ρ`
`F` and `ρ` must be of the same size.
"""
function eig_cf(N, dx, ρ, F, k::Number, bcs::Tuple{T,U}; η=0,
                ncf::Int=1, ka::Number=0, kb::Number=0,
                coordinate_system::Symbol=:cart, h=1) where U<:BoundaryCondition where T<:BoundaryCondition

    ∇², S = laplacian(N, dx, bcs; ka=ka, kb=kb, coordinate_system=coordinate_system, h=h)

    ɛk² = spdiagm(0 => ɛ[:]*k^2)
    Fk² = spdiagm(0 => -F[:]*k^2)

    decomp, history = partialschur(shift_and_invert(∇²+ɛk², Fk², η), nev=ncf)
    history.converged ? nothing : @warn "incomplete convergence: only $(h.nconverged) of $(h.nev) evecs converged"
    η, u = η .+ 1 ./decomp.eigenvalues, decomp.Q
    for i ∈ 1:ncf
        𝒩² = sum(u[:,i].*F[:].*conj(u[:,i]))*prod(dx)
        @views u[:,i] = u[:,i]/sqrt(𝒩²)
    end
    return η, u
end
function eig_cf(dx, ρ, F, k::Number, bcs::Tuple{T,U}; η=0,
                ncf::Int=1, ka::Number=0, kb::Number=0,
                coordinate_system::Symbol=:cart, h=ones(N...)) where U<:BoundaryCondition where T<:BoundaryCondition
    return eig_cf([size(ρ)...], dx, ρ, F, k, bcs; η=η, ncf=ncf, ka=ka, kb=kb, coordinate_system=coordinate_system, h=h)
end


end # module
