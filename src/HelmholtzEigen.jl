module HelmholtzEigen

include("ArnoldiWrapper.jl"); using .ArnoldiWrapper
using ..CoordinateSystems,
..BoundaryConditions,
..DifferentialOperators,
..SimulationDefinition,
ArnoldiMethod,
LinearAlgebra,
SparseArrays,
NonlinearEigenproblems,
Distributed,
LinearAlgebra,
Random,
ProgressMeter

export eig_kl,
eig_cf,
eig_knl

function eig_kl(args...; disp_opt::Bool=false, kwargs...)
    decomp, history = partialschur(args...; cf=false, kwargs...)
    disp_opt ? println(history) : nothing
    @assert history.converged history
    return decomp.eigenvalues, Array(decomp.Q)
end
function eig_cf(args...; disp_opt::Bool=false, kwargs...)
    decomp, history = partialschur(args...; cf=true, kwargs...)
    disp_opt ? println(history) : nothing
    @assert history.converged history
    return decomp.eigenvalues, Array(decomp.Q)
end
function eig_knl(args...; quad_n::Int=100, disp_opt::Bool=false, method::Symbol=:contour_beyn, nk::Int=3, kwargs...)
    if method==:contour_beyn
        displaylevel = disp_opt ? 1 : 0
        k, ψ = contour_beyn(args...; N=quad_n, displaylevel=displaylevel, k=nk, kwargs...)
    else
        throw(ArgumentError("unrecognized method $method"))
    end
    return k, ψ
end


function ArnoldiMethod.partialschur(sim::Simulation, k::Number, ka::Number=0, kb::Number=0, lattice_index::Int=0; cf::Bool=false, η::Number=0, kwargs...)
    N,args = sim.dis.N,(sim.dis.coordinate_system,sim.bnd.∂Ω,sim.bnd.bc,sim.bnd.bl)
    (∇²,fs,s),defs = laplacian(N,args...,k,ka,kb,lattice_index)
    κ = lattice_index==1 ? ka : kb
    A = spzeros(eltype(∇²[1]),prod(N),prod(N))
    for i ∈ eachindex(∇²)
        A += ∇²[i]*fs[i](k,κ)
    end
    A = (A + SparseMatrixCSC(transpose(A)))/2
    B = s*spdiagm(0=>-sim.sys.ε[:])
    if !cf
        σ=k^2
    else
        A -= B*k^2
        B = s*spdiagm(0=>-sim.sys.F[:]*k^2)
        σ=η
    end

    decomp, history = partialschur(A, B, σ; diag_inv_B=true, kwargs...)

    !cf ? decomp.eigenvalues[:] = sqrt.(decomp.eigenvalues[:]) : nothing

    # Normalize wavefunctions according to (ψ₁,ψ₂)=δ₁₂, which requires transformed ε or F
    normalize!(sim,decomp.Q,B)
    return decomp, history
end


function NonlinearEigenproblems.contour_beyn(sim::Simulation, k::Number, ka::Number=0, kb::Number=0, lattice_index::Int=0; quad_n::Int=100, displaylevel::Int=0, kwargs...)
    nep = SPMF_NEP(sim, k; check_consistency=false)
    if displaylevel>0
        k, ψ = contour_beyn(nep, true; N=quad_n, σ=k, kwargs...)
    else
        k, ψ = contour_beyn(nep; N=quad_n, σ=k, displaylevel=0, kwargs...)
    end
    normalize!(sim,ψ,nep.A[end])
    return k, ψ
end

function NonlinearEigenproblems.SPMF_NEP(sim::Simulation, k::Number, ka::Number=0, kb::Number=0, lattice_index::Int=0; kwargs...)
    N,brgs = sim.dis.N,(sim.dis.coordinate_system,sim.bnd.∂Ω,sim.bnd.bc,sim.bnd.bl)
    (∇²,fs,s),defs = laplacian(N,brgs...,k,ka,kb,lattice_index)
    for i ∈ eachindex(∇²)
        ∇²[i] = (∇²[i] + SparseMatrixCSC(transpose(∇²[i])))/2
    end
    B = s*spdiagm(0=>sim.sys.ε[:])
    push!(∇²,B)
    push!(fs,k->k^2)
    return SPMF_NEP(∇²,fs; kwargs...)
end

function normalize!(sim::Simulation,ψ,B)
    dx = sim.dis.dx
    for i ∈ 1:size(ψ,2)
        𝒩² = sum((ψ[:,i].^2).*diag(B))*(isinf(dx[1]) ? 1 : dx[1])*(isinf(dx[2]) ? 1 : dx[2])
        ψ[:,i] /= sqrt(𝒩²)*exp(complex(0,angle(ψ[end÷2-1,i])))
    end
    return nothing
end


NonlinearEigenproblems.contour_beyn(nep::NEP,disp_opt::Bool;params...)=contour_beyn(ComplexF64,nep,disp_opt;params...)
function NonlinearEigenproblems.contour_beyn(::Type{T},
                         nep::NEP,
                         disp_opt::Bool;
                         tol::Real=eps(real(T))*100,
                         rank_tol::Real=sqrt(eps()),
                         σ::Number=zero(complex(T)),
                         displaylevel::Integer=0,
                         linsolvercreator::Function=backslash_linsolvercreator,
                         k::Integer=3, # Number of eigenvals to compute
                         radius=1, # integration radius
                         parallel::Bool = nprocs()>1,
                         quad_method::Symbol= parallel ? :ptrapz_parallel : :ptrapz, # which method to run. :quadg, :quadg_parallel, :quadgk, :ptrapz, :ptrapz_parallel
                         N::Integer=1000  # Nof quadrature nodes
                         )where{T<:Number}

    n=size(nep,1);
    Random.seed!(10); # Reproducability
    Vh=Array{T,2}(randn(real(T),n,k)) # randn only works for real
    temp = Array{T,2}(undef,n,k)

    length(radius)==1 ? radius=(radius,radius) : nothing
    function ggp(t)
        sc = sincos(t)
        g = complex(radius[1]*sc[1],radius[2]*sc[2])
        gp = complex(-radius[1]*sc[2],radius[2]*sc[1])
        return g,gp
    end

    if (k>n)
        println("k=",k," n=",n);
        error("Cannot compute more eigenvalues than the size of the NEP with contour_beyn()");
    end

    function local_linsolve(λ::TT,V::Matrix{TT}) where {TT<:Number}
        @ifd(print("."))
        local M0inv::LinSolver = linsolvercreator(nep,λ+σ);
        # This requires that lin_solve can handle rectangular
        # matrices as the RHS
        return lin_solve(M0inv,V);
    end

    # Constructing integrands
    function Tv!(λ,t)
        t[:] = local_linsolve(T(λ),Vh)
        return nothing
    end
    function f!(λ,dλ,t)
        Tv!(λ,t)
        t[:] = t*dλ
        return nothing
    end
    @ifd(print("Computing integrals"))

    local A0,A1
    pg = Progress(N; dt=.1, desc="Contour integration...")
    if (quad_method == :quadg_parallel)
        println(" using quadg_parallel")
        #A0=quadg_parallel(f1,0,2*pi,N);
        #A1=quadg_parallel(f2,0,2*pi,N);
        error("disabled");
    elseif (quad_method == :quadg)
        println(" using quadg")
        #A0=quadg(f1,0,2*pi,N);
        #A1=quadg(f2,0,2*pi,N);
        error("disabled");
    elseif (quad_method == :ptrapz)
        println(" using ptrapz")
        (A0,A1)=ptrapz((ggp,f!),0,2*pi,N,temp,pg);
    elseif (quad_method == :ptrapz_parallel)
        println(" using ptrapz_parallel")
        channel = RemoteChannel(()->Channel{Bool}(pg.n+1),1)
        @sync begin
            @async begin
                (A0,A1)=ptrapz_parallel((ggp,f!),0,2*pi,N,temp,channel);
                put!(channel,false)
            end
            @async while take!(channel)
                next!(pg)
            end
        end
    elseif (quad_method == :quadgk)
        println(" using quadgk")
        # A0,tmp=quadgk(f1,0,2*pi,pg...,1,reltol=tol);
        # A1,tmp=quadgk(f2,0,2*pi,pg...,2,reltol=tol);
        error("disabled");
    else
        println()
        error("Unknown quadrature method: "*String(quad_method));
    end
    # Don't forget scaling
    A0[:,:] = A0 ./(2im*pi);
    A1[:,:] = A1 ./(2im*pi);

    println("Computing SVD prepare for eigenvalue extraction ")
    V,S,W = svd(A0)
    max_ind = findlast(S./minimum(S).>1/rank_tol)
    if isnothing(max_ind)
        @warn "Rank not revealed, contour likely encloses too many eigenvalues. Try increasing k"
        max_ind=k
    end
    V0, W0, S0 = V[:,1:max_ind], W[:,1:max_ind], S[1:max_ind]
    B = (V0'*A1*W0)*diagm(0 => 1 ./S0[:])

    println("Computing eigenvalues ")
    λ,v=eigen(B)
    λ[:] = λ .+ σ

    println("Computing eigenvectors ")
    return (λ,V0*v)
end

# Trapezoidal rule for a periodic function f
function ptrapz((ggp,f!),a,b,N,temp,pg)
    h = (b-a)/N
    t = range(a, stop = b-h, length = N)
    S = [zero(temp),zero(temp)]
    for i = 1:N
        g, gp = ggp(t[i])
        f!(g,gp,temp)
        S[1][:] = S[1] + temp
        S[2][:] = S[2] + temp*g
        next!(pg)
    end
    return h*S[1], h*S[2];
end
# parallel Trapezoidal rule for a periodic function f
function ptrapz_parallel((ggp,f!),a,b,N,temp,channel)
    h = (b-a)/N
    t = range(a, stop = b-h, length = N)
    S = [zero(temp),zero(temp)]
    S = @distributed (+) for i = 1:N
        temp = copy(temp)
        g, gp = ggp(t[i])
        f!(g,gp,temp)
        put!(channel,true)
        [temp,temp*g]
    end
    return h*S[1], h*S[2];
end

end # module
