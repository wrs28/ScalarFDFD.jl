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

number of lattice sites `N`, lattice spacing `dx`, density function `ρ`, boundary conditions `bcs` (see `DifferentialOperators.ezbc`),
number of frequencies `nk`, Floquet/Bloch wavenumbers `ka`,`kb`.
"""
function eig_kl(N::Array{Int,1},
                dx::Array{Float64,1},
                ρ::Array{Tρ,2},
                k::Number,
                bcs::Tuple{Tbc1,Tbc2};
                nk::Int=1,
                ka::Number=0,
                kb::Number=0,
                coordinate_system::Symbol=:cart,
                h::Array{Th,2}=ones(size(ρ,1)+1, size(ρ,2)+1)
                ) where Tρ<:Number where Tbc1<:BoundaryCondition where Tbc2<:BoundaryCondition where Th<:Number

    ∇², S = laplacian(N, dx, bcs; ka=ka, kb=kb, coordinate_system=coordinate_system, h=h)
    decomp, history = partialschur(shift_and_invert(∇², spdiagm(0=>-ρ[:]), k^2, diag_inv_B=true), nev=nk)
    history.converged ? nothing : @warn "incomplete convergence: only $(h.nconverged) of $(h.nev) evecs converged"
    k², ψ = k^2 .+ 1 ./decomp.eigenvalues, decomp.Q
    for i ∈ 1:nk
        𝒩² = sum(abs2.(ψ[:,i].*ρ[:]))*(isinf(dx[1]) ? 1 : dx[1])*(isinf(dx[2]) ? 1 : dx[2])
        @views ψ[:,i] = ψ[:,i]/sqrt(𝒩²)/exp(complex(0,angle(ψ[end÷2,i])))
    end
    return sqrt.(k²), convert(Array,ψ)
end


"""
    η, u =  eig_cf(N, dx, ρ, k, bcs; nk=1, ka=0, kb=0, coordinate_system=:cart, h=1)

`F` and `ρ` must be of the same size.
"""
function eig_cf(N::Array{Int,2},
                dx::Array{Float64,1},
                ρ::Array{Tρ,2},
                F::Array{Tf,2},
                k::Number,
                bcs::Tuple{Tbc1,Tbc2};
                η::Number=0,
                ncf::Int=1,
                ka::Number=0,
                kb::Number=0,
                coordinate_system::Symbol=:cart,
                h::Array{Th,2}=ones(size(ρ,1)+1, size(ρ,2)+1)
                ) where Tρ<:Number where Tf<:Number where Tbc1<:BoundaryCondition where Tbc2<:BoundaryCondition where Th<:Number

    ∇², S = laplacian(N, dx, bcs; ka=ka, kb=kb, coordinate_system=coordinate_system, h=h)
    ɛk² = spdiagm(0 => ɛ[:]*k^2)
    Fk² = spdiagm(0 => -F[:]*k^2)
    decomp, history = partialschur(shift_and_invert(∇²+ɛk², Fk², η), nev=ncf)
    history.converged ? nothing : @warn "incomplete convergence: only $(h.nconverged) of $(h.nev) evecs converged"
    η, u = η .+ 1 ./decomp.eigenvalues, decomp.Q
    for i ∈ 1:ncf
        𝒩² = sum(u[:,i].*F[:].*conj(u[:,i]))*(isinf(dx[1]) ? 1 : dx[1])*(isinf(dx[2]) ? 1 : dx[2])
        @views u[:,i] = u[:,i]/sqrt(𝒩²)
    end
    return η,  convert(Array,u)
end


end # module
