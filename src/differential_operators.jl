# TODO: open boundaries in laplacians_with_bc

function ⊗(A,B)
    kron(A,B)
end

################################################################################
####### GRADIENTS
################################################################################
"""
    ∇₁, ∇₂ =  grad(sim; symmetric=false)

2-dim gradients with `sim.dis.N` points, lattice spacing `sim.dis.dx`.
"""
function grad(sim::Simulation; symmetric=false)
    return grad(sim.dis; symmetric=symmetric)
end


"""
    ∇₁, ∇₂ =  grad(dis; symmetric=false)

2-dim gradients with `dis.N` points, lattice spacing `dis.dx`.
"""
function grad(dis::Discretization; symmetric=false)
    return grad(dis.N, dis.dx; symmetric=symmetric)
end


"""
    ∇ =  grad(N::Int, dx::Float64; symmetric=false)

`symmetric=false`: 1-dim forward/backward Gradient (depending on interpretation of ∇ψ) with `N` points,
lattice spacing `dx`.

`symmetric=true`: 1-dim symmetric Gradient (depending on interpretation of ∇ψ) with `N` points,
lattice spacing `dx`.
"""
function grad(N::Int, dx; symmetric=false)
    if symmetric
        I₁ = Array(2:N)
        J₁ = Array(1:N-1)
        V₁ = fill(ComplexF64(-1/2dx), N-1)

        I₂ = Array(1:N-1)
        J₂ = Array(2:N)
        V₂ = fill(ComplexF64(+1/2dx), N-1)

        ∇ = sparse(vcat(I₁,I₂), vcat(J₁,J₂), vcat(V₁,V₂), N, N, +)
    else
        I₁ = Array(1:N-1)
        J₁ = Array(1:N-1)
        V₁ = fill(ComplexF64(-1/dx), N-1)

        I₂ = Array(1:N-1)
        J₂ = Array(2:N)
        V₂ = fill(ComplexF64(+1/dx), N-1)

        ∇ = sparse(vcat(I₁,I₂), vcat(J₁,J₂), vcat(V₁,V₂), N-1, N, +)
    end
    return ∇
end



"""
    ∇₁, ∇₂ =  grad(N::Array{Int}, dx::Array{Float64}; symmetric=false)

2-dim gradients with `N[1]`, `N[2]` points, lattice spacing `dx[1], dx[2]`.
"""
function grad(N::Array{Int}, dx; symmetric=false)
    ∇1 = grad(N[1],dx[1]; symmetric=symmetric)
    ∇2 = grad(N[2],dx[2]; symmetric=symmetric)

    𝕀1 = sparse(complex(1.,0)I, N[1], N[1])
    𝕀2 = sparse(complex(1.,0)I, N[2], N[2])

    ∇₁ = 𝕀2 ⊗ ∇1
    ∇₂ = ∇2 ⊗ 𝕀1

    return ∇₁, ∇₂
end


################################################################################
####### LAPLACIANS
################################################################################
"""
    ∇² = laplacian(sim, k; ka=0, kb=0)
"""
function laplacian(sim::Simulation, k; ka=0, kb=0)

    ∇₁², ∇₂² = ScalarFDFD.laplacians_sans_bc(sim, k)
    ScalarFDFD.laplacians_with_bc!(∇₁², ∇₂², sim)

    C1 = sim.bnd.weights[1]
    C1ᵀ = sim.bnd.weights[2]
    N1a = sim.bnd.weights[3]
    N1b = sim.bnd.weights[4]
    N1aᵀ = sim.bnd.weights[5]
    N1bᵀ = sim.bnd.weights[6]
    C2 = sim.bnd.weights[7]
    C2ᵀ = sim.bnd.weights[8]
    N2a = sim.bnd.weights[9]
    N2b = sim.bnd.weights[10]
    N2aᵀ = sim.bnd.weights[11]
    N2bᵀ = sim.bnd.weights[12]

    𝕀1 = sparse(complex(1.,0)I, sim.dis.N[1], sim.dis.N[1])
    𝕀2 = sparse(complex(1.,0)I, sim.dis.N[2], sim.dis.N[2])

    ϕ1 = ka*sim.lat.a*N1a + kb*sim.lat.b*N1b
    ϕ1ᵀ = ka*sim.lat.a*N1aᵀ + kb*sim.lat.b*N1bᵀ
    ϕ2 = ka*sim.lat.a*N2a + kb*sim.lat.b*N2b
    ϕ2ᵀ = ka*sim.lat.a*N2aᵀ + kb*sim.lat.b*N2bᵀ

    C1 = C1.*exp.(+1im*ϕ1) + C1ᵀ.*exp.(-1im*ϕ1ᵀ)
    C2 = C2.*exp.(+1im*ϕ2) + C2ᵀ.*exp.(-1im*ϕ2ᵀ)

    # ϕ1 = ka*sim.lat.a*N1a.nzval + kb*sim.lat.b*N1b.nzval
    # ϕ1ᵀ = ka*sim.lat.a*N1aᵀ.nzval + kb*sim.lat.b*N1bᵀ.nzval
    # ϕ2 = ka*sim.lat.a*N2a.nzval + kb*sim.lat.b*N2b.nzval
    # ϕ2ᵀ = ka*sim.lat.a*N2aᵀ.nzval + kb*sim.lat.b*N2bᵀ.nzval

    # i1, j1, k1 = findnz(Cₐ)
    # i2, j2, k2 = findnz(Cₐᵀ)
    # C1 = sparse(i1,j1, k1.*exp.(+1im*ϕ1), prod(sim.dis.N), prod(sim.dis.N)) +
         # sparse(i2,j2, k2.*exp.(-1im*ϕ1ᵀ), prod(sim.dis.N), prod(sim.dis.N))

    # i1, j1, k1 = findnz(Cᵦ)
    # i2, j2, k2 = findnz(Cᵦᵀ)
    # C2 = sparse(i1,j1, k1.*exp.(+1im*ϕ2), prod(sim.dis.N), prod(sim.dis.N)) +
         # sparse(i2,j2, k2.*exp.(-1im*ϕ2ᵀ), prod(sim.dis.N), prod(sim.dis.N))

    return (
    (𝕀2 ⊗ ∇₁²) + (∇₂² ⊗ 𝕀1) + C1 + C2
    )
end


"""
    ∇₁², ∇₂² = laplacians_sans_bc(sim, k)
"""
function laplacians_sans_bc(sim::Simulation, k)
    N = sim.dis.N
    dx = sim.dis.dx
    Σ = σ(sim)
    S1 = [sparse(1:N[j]-1, 1:N[j]-1, Vector{ComplexF64}(undef,N[j]-1), N[j]-1, N[j]-1) for i ∈ 1:2, j ∈ 1:2]
    S2 = [sparse(1:N[j], 1:N[j], Vector{ComplexF64}(undef,N[j]), N[j], N[j]) for i ∈ 1:2, j ∈ 1:2]
    SA = [sparse(1:N[j], 1:N[j], Vector{ComplexF64}(undef,N[j]), N[j], N[j]) for i ∈ 1:2, j ∈ 1:2]
    for r ∈ CartesianIndices(S1)
        j = r[2]
        if sim.bnd.bl[r] ∈ [:pml_out, :pml_in]
            S1[r] = sparse(1:N[j]-1, 1:N[j]-1, 1 ./(1 .+ 1im*(Σ[r][1:end-1] + Σ[r][2:end])/real(2k)), N[j]-1, N[j]-1)
            S2[r] = sparse(1:N[j], 1:N[j], 1 ./(1 .+ 1im*Σ[r]/real(k)), N[j], N[j])
            SA[r] = sparse(complex(1.,0)I, N[j], N[j])
        else
            S1[r] = sparse(complex(1.,0)I, N[j]-1, N[j]-1)
            S2[r] = sparse(complex(1.,0)I, N[j], N[j])
            SA[r] = sparse(1:N[j], 1:N[j], 1 ./(1 .+ 1im*Σ[r]/real(k)), N[j], N[j])
        end
    end

    ∇₁ = grad(N[1], dx[1])
    ∇₂ = grad(N[2], dx[2])

    ∇₁² = -(SA[1,1]*SA[1,1]*SA[2,1]*SA[2,1]*S2[1,1]*S2[2,1]*transpose(∇₁)*S1[1,1]*S1[2,1]*∇₁)
    ∇₂² = -(SA[1,2]*SA[1,2]*SA[2,2]*SA[2,2]*S2[1,2]*S2[2,2]*transpose(∇₂)*S1[1,2]*S1[2,2]*∇₂)

    return ∇₁², ∇₂²
end


"""
    laplacians_with_bc!(∇₁², ∇₂², sim)
"""
function laplacians_with_bc!(∇₁², ∇₂², sim::Simulation)
    bc = sim.bnd.bc
    periodic_boundary_weights!(sim)

    Nx = sim.dis.N[1]
    dx = sim.dis.dx[1]
    dx² = dx^2
    ind = [1, Nx]
    if Nx > 1
        for i ∈ 1:2
            if bc[i,1] == :d
                ∇₁²[ind[i],ind[i]] += -2/dx²
            elseif bc[i,1] == :n
                ∇₁²[ind[i],ind[i]] += 0
            elseif bc[i,1] == :p
                ∇₁²[ind[i],ind[i]] += -1/dx²
            end
        end
    end

    Ny = sim.dis.N[2]
    dy = sim.dis.dx[2]
    dy² = dy^2
    ind = [1, Ny]
    if Ny > 1
        for i ∈ 1:2
            if bc[i,2] == :d
                ∇₂²[ind[i],ind[i]] += -2/dy²
            elseif bc[i,2] == :n
                ∇₂²[ind[i],ind[i]] += 0
            elseif bc[i,2] == :p
                ∇₂²[ind[i],ind[i]] += -1/dy²
            end
        end
    end
    return nothing
end


################################################################################
####### AUXILLIARIES
################################################################################
"""
    periodic_boundary_weights!(sim)

compute periodic boundary weights, saves to sim.bnd.weights
"""
function periodic_boundary_weights!(sim::Simulation)

    try
        sim.bnd.weights[1]
    catch
        bc = sim.bnd.bc
        dx = sim.dis.dx[1]; dx² = dx^2
        dy = sim.dis.dx[2]; dy² = dy^2

        Cₐ = spzeros(Float64,prod(sim.dis.N),prod(sim.dis.N))
        Cₐᵀ = spzeros(Float64,prod(sim.dis.N),prod(sim.dis.N))
        Cᵦ = spzeros(Float64,prod(sim.dis.N),prod(sim.dis.N))
        Cᵦᵀ = spzeros(Float64,prod(sim.dis.N),prod(sim.dis.N))
        N1a = spzeros(Float64,prod(sim.dis.N),prod(sim.dis.N))
        N1b = spzeros(Float64,prod(sim.dis.N),prod(sim.dis.N))
        N2a = spzeros(Float64,prod(sim.dis.N),prod(sim.dis.N))
        N2b = spzeros(Float64,prod(sim.dis.N),prod(sim.dis.N))

        if :p ∈ bc[:,1] && bc[1,1] == bc[2,1]
            I, J, V, Na, Nb =  periodic_boundary_weights(sim, 1)
            Cₐ   = sparse(I, J, V/dx², prod(sim.dis.N), prod(sim.dis.N))
            Cₐᵀ  = sparse(J, I, V/dx², prod(sim.dis.N), prod(sim.dis.N))
            N1a  = sparse(I, J, Na, prod(sim.dis.N), prod(sim.dis.N), mean)
            N1aᵀ = sparse(J, I, Na, prod(sim.dis.N), prod(sim.dis.N), mean)
            N1b  = sparse(I, J, Nb, prod(sim.dis.N), prod(sim.dis.N), mean)
            N1bᵀ = sparse(J, I, Nb, prod(sim.dis.N), prod(sim.dis.N), mean)
        elseif :p ∈ bc[:,1]
            throw(ArgumentError("only one boundary of dimension 1 is periodic, must be both or none"))
        end

        if :p ∈ bc[:,2] && bc[1,2] == bc[2,2]
            I, J, V, Na, Nb = periodic_boundary_weights(sim, 2)
            Cᵦ  = sparse(I, J, V/dy², prod(sim.dis.N), prod(sim.dis.N))
            Cᵦᵀ = sparse(J, I, V/dy², prod(sim.dis.N), prod(sim.dis.N))
            N2a  = sparse(I, J, Na, prod(sim.dis.N), prod(sim.dis.N), mean)
            N2aᵀ = sparse(J, I, Na, prod(sim.dis.N), prod(sim.dis.N), mean)
            N2b  = sparse(I, J, Nb, prod(sim.dis.N), prod(sim.dis.N), mean)
            N2bᵀ = sparse(J, I, Nb, prod(sim.dis.N), prod(sim.dis.N), mean)
        elseif :p ∈ bc[:,2]
            throw(ArgumentError("only one boundary of dimension 2 is periodic, must be both or none"))
        end

        sim.bnd.weights[1] = Cₐ
        sim.bnd.weights[2] = Cₐᵀ
        sim.bnd.weights[3] = N1a
        sim.bnd.weights[4] = N1b
        sim.bnd.weights[5] = N1aᵀ
        sim.bnd.weights[6] = N1bᵀ
        sim.bnd.weights[7] = Cᵦ
        sim.bnd.weights[8] = Cᵦᵀ
        sim.bnd.weights[9] = N2a
        sim.bnd.weights[10] = N2b
        sim.bnd.weights[11] = N2aᵀ
        sim.bnd.weights[12] = N2bᵀ
    end
    return nothing
end

"""

"""
function remove_boundary_weights!(sim::Simulation)
    for i ∈ 1:4
        pop!(sim.bnd.weights)
    end
    return nothing
end


"""
    I, J, V = periodic_boundary_weights(sim, dim)
"""
function periodic_boundary_weights(sim::Simulation, dim)
    if !(iszero(sim.lat.α) || iszero(sim.lat.β-π/2))
        throw(ArgumentError("this scheme does not work for angled periodic lattices.
    This lattice has α=$(sim.lat.α) and β=$(sim.lat.β)"))
    end

    sim = deepcopy(sim)
    lattice = Bravais(sim.lat; :x0=>0, :y0=>0)
    ∂Ω = sim.bnd.∂Ω
    ∂Ω_tr = sim.bnd.∂Ω_tr
    bc = sim.bnd.bc
    bl = sim.bnd.bl
    bl_depth = sim.bnd.bl_depth

    a = lattice.a
    b = lattice.b
    α = lattice.α
    β = lattice.β

    N = sim.dis.N

    if dim == 1
        if !isinf(a)
            ∂Ω[1,1] = min(0, a*(cos(α)-sin(α)*cot(β)))
            ∂Ω[2,1] = max(0, a*(cos(α)-sin(α)*cot(β)))
        else
            return Array{Int}(undef,0), Array{Int}(undef,0), Array{Float64}(undef,0), Array{Float64}(undef,0)
        end
        ∂Ω[2,2] = ∂Ω[2,2]-∂Ω[1,2]
        ∂Ω[1,2] = 0
    elseif dim == 2
        ∂Ω[2,1] = ∂Ω[2,1]-∂Ω[1,1]
        ∂Ω[1,1] = 0
        if !isinf(b)
            ∂Ω[1,2] = min(0, b*sin(β))
            ∂Ω[2,2] = max(0, b*sin(β))
        else
            return Array{Int}(undef,0), Array{Int}(undef,0), Array{Float64}(undef,0), Array{Float64}(undef,0)
        end
    else
        throw(ArgumentError("invalid dimensions $(dim)"))
    end

    bnd = Boundary(∂Ω, bc, bl, bl_depth)
    sim = Simulation(sim; :bnd => bnd, :lat => lattice )
    x = sim.dis.x[1]; y = sim.dis.x[2]; dx = sim.dis.dx[1]; dy = sim.dis.dx[2]
    lattice = sim.lat; N = sim.dis.N

    if dim == 1
        P = ScalarFDFD.bravais_coordinates.(x[1]-dx, y, Ref(lattice))
    else
        P = ScalarFDFD.bravais_coordinates.(x, y[1]-dy, Ref(lattice))
    end

    p1 = Array{Float64}(undef, N[mod1(dim+1,2)])
    p2 = Array{Float64}(undef, N[mod1(dim+1,2)])
    for i ∈ eachindex(P)
        p1[i] = P[i][1]
        p2[i] = P[i][2]
    end
    Ma = -floor.(Int, p1/sim.lat.a)
    Mb = -floor.(Int, p2/sim.lat.b)
    X = sim.lat.v1[1]*(p1 + Ma*sim.lat.a) + sim.lat.v2[1]*(p2 + Mb*sim.lat.b)
    Y = sim.lat.v1[2]*(p1 + Ma*sim.lat.a) + sim.lat.v2[2]*(p2 + Mb*sim.lat.b)

    Ma += -floor.(Int, X/(sim.bnd.∂Ω[2,1]-sim.bnd.∂Ω[1,1]))
    Mb += -floor.(Int, Y/(sim.bnd.∂Ω[2,2]-sim.bnd.∂Ω[1,2]))

    x_inds1 = floor.(Int, X/dx .+ 1/2)
    x_inds2 = x_inds1 .+ 1

    y_inds1 = floor.(Int, Y/dy .+ 1/2)
    y_inds2 = y_inds1 .+ 1

    Cx1 = abs.(X/dx .+ 1/2 - x_inds2)
    Cx2 = abs.(X/dx .+ 1/2 - x_inds1)
    cx1 = Cx1./(Cx1+Cx2)
    cx2 = Cx2./(Cx1+Cx2)

    Cy1 = abs.(Y/dy .+ 1/2 - y_inds2)
    Cy2 = abs.(Y/dy .+ 1/2 - y_inds1)
    cy1 = Cy1./(Cy1+Cy2)
    cy2 = Cy2./(Cy1+Cy2)

    q = Array{Int}(undef,2)
    r = Array{Int}(undef,2)
    s = Array{Float64}(undef,2)

    t = Array{Int}(undef,2)
    u = Array{Int}(undef,2)
    v = Array{Float64}(undef,2)

    j = Array{Int}(undef,4)
    k = Array{Int}(undef,4)
    l = Array{Float64}(undef,4)

    I = Array{Int}(undef, 4N[mod1(dim+1,2)])
    J = Array{Int}(undef, 4N[mod1(dim+1,2)])
    V = Array{Float64}(undef, 4N[mod1(dim+1,2)])
    Na = Array{Float64}(undef, 4N[mod1(dim+1,2)])
    Nb = Array{Float64}(undef, 4N[mod1(dim+1,2)])

    # pg = Progress(N[mod1(dim+1,2)], PROGRESS_UPDATE_TIME::Float64, "periodic boundaries ")
    for i ∈ 1:N[mod1(dim+1,2)]

        if dim == 1
            ind_x = 1
            ind_y = i
        else
            ind_x = i
            ind_y = 1
        end

        q[1] = ind_x
        q[2] = ind_x

        r[1] = mod1(x_inds1[i],N[1])
        r[2] = mod1(x_inds2[i],N[1])

        s[1] = cx1[i]
        s[2] = cx2[i]

        t[1] = ind_y
        t[2] = ind_y

        u[1] = mod1(y_inds1[i],N[2])
        u[2] = mod1(y_inds2[i],N[2])

        v[1] = cy1[i]
        v[2] = cy2[i]


        if 1 ∈ N
            j[1:2], k[1:2], l[1:2] = findnz( sparse(t, u, v, N[2], N[2]) ⊗ sparse(q, r, s, N[1], N[1]) )
            j[3:4] .= k[3:4] .= 1
            l[3:4] .= 0
        else
            j[:], k[:], l[:] = findnz( sparse(t, u, v, N[2], N[2]) ⊗ sparse(q, r, s, N[1], N[1]) )
        end

        I[(4(i-1)+1):(4(i-1)+4)] = j
        J[(4(i-1)+1):(4(i-1)+4)] = k
        V[(4(i-1)+1):(4(i-1)+4)] = l
        Na[(4(i-1)+1):(4(i-1)+4)] .= Ma[i]
        Nb[(4(i-1)+1):(4(i-1)+4)] .= Mb[i]

        # next!(pg)
    end
    return I, J, V, Na, Nb
end


"""
    s₁, s₂ = σ(sim)

conductivity for absorbing layer (PML or not) in dimensions 1 and 2.
"""
function σ(sim::Simulation)

    α = Array{ComplexF64}(undef,2,2)
    for i ∈ CartesianIndices(α)
        if sim.bnd.bl[i] !== :none
            α[i] = -(1/4)*(POWER_LAW+1)*exp(complex(0,SCALING_ANGLE))*log(EXTINCTION::Float64)/(sim.bnd.bl_depth[i]^(float(POWER_LAW)+1))
        else
            α[i] = 0
        end
        if sim.bnd.bl[i] ∈ [:abs_in, :pml_in]
            α[i] = flipsign(conj(α[i]),-1)
        end
    end

    Σ = Array{Array{ComplexF64,1},2}(undef,2,2)
    for r ∈ CartesianIndices(Σ)
        i = r[1]
        j = r[2]
        Σ[i,j] = zeros(ComplexF64,length(sim.dis.x[j]))
        if sim.bnd.bl[r] !== :none
            for k ∈ eachindex(sim.dis.x[j])
                if sign(sim.bnd.∂Ω_tr[i,j] - sim.dis.x[j][k])*(-1)^i ≤ 0
                     Σ[i,j][k] = α[i,j]*(abs(sim.dis.x[j][k]-sim.bnd.∂Ω_tr[i,j]))^float(POWER_LAW)
                end
            end
        end
    end

    return Σ
end
