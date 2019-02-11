module ArnoldiFormat

export shift_and_invert

using ArnoldiMethod,
LinearAlgebra,
LinearMaps,
SparseArrays


"""
    type ShiftAndInvert{T,U,V}

Fieds:
`A::T`, `B::U`, `temp::V`.

Used to contain arrays for function `shift_and_invert` for use in `ArnoldiMethod.partialschur`
"""
struct ShiftAndInvert{T,U,V}
    A_lu::T
    B::U
    temp::V
end


"""
    ShiftAndInvert(y,x)

A\\B*x = y, where `A=ShiftAndInvert.M.A_lu`, `B=ShiftAndInvert.M.B`, `x=ShiftAndInvert.M.x`
"""
function (M::ShiftAndInvert)(y,x)
    mul!(M.temp, M.B, x)
    ldiv!(y, M.A_lu, M.temp)
    return nothing
end


"""
    a = shift_and_invert(A, B, σ; diag_inv_B=false)

create a LinearMap object to feed to ArnoldiMethod.partialschur which transforms `Ax=λBx` into `(A-σB)⁻¹Bx=x/(λ-σ)`.

Set `diag_inv_B`=true if `B` is both diagonal and invertible, then it is easy to compute `B⁻¹`, so instead returns linear map `(B⁻¹A-σI)⁻¹`, which has same evals as above.

`A` and `B` must both be sparse or both dense. `A`, `B`, `σ` need not have common element type.
"""
function shift_and_invert(A::S, B::T, σ::Number; diag_inv_B::Bool=true, issymmetric::Bool=false) where {S,T}

    onetype = one(eltype(A))*one(eltype(B))*one(σ)
    type = typeof(onetype)
    A, B = onetype*A, onetype*B
    @assert supertype(typeof(A))==supertype(typeof(B)) "typeof(A)=$(typeof(A)), typeof(B)=$(typeof(B)). Either both sparse or both dense"
    if diag_inv_B
        if T<:AbstractSparseArray
            diag_fun = spdiagm
            matrix_fun = sparse
        else
            diag_fun = Diagonal
            matrix_fun = Matrix
        end
        α = diag_fun(0=>map(x->1/x,diag(B)))*A-σ*I
        β = matrix_fun(onetype*I,size(B))
    else
        α = A-σ*B
        β = onetype*B
    end
    a = ShiftAndInvert(lu(α), β, Vector{type}(undef, size(A,1)))

    return LinearMap{type}(a, size(A,1); ismutating=true, issymmetric=issymmetric)
end


"""
    a = shift_and_invert(A, σ; diag_inv_B=false)
"""
function shift_and_invert(A::S, σ::Number; kwargs...) where S
    onetype = one(eltype(A))*one(σ) ; type = typeof(onetype)
    if S<:AbstractSparseArray
        return shift_and_invert(A,sparse(onetype*I,size(A)...), σ; diag_inv_B=true, kwargs...)
    else
        return shift_and_invert(A,Matrix(onetype*I,size(A)...), σ; diag_inv_B=true, kwargs...)
    end
end

function ArnoldiFormat.partialeigen(decomp,σ)
    λ, ψ = partialeigen(decomp)
    return (σ.+ 1 ./λ), ψ
end

end # module ArnoldiFormat
