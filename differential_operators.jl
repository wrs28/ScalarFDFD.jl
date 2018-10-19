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
    ∇₁², ∇₂² = laplacians_sans_bc(sim, k)
    laplacians_with_bc!(∇₁², ∇₂², sim)
    C₁  = sim.bnd.weights[1]
    C₁ᵀ = sim.bnd.weights[2]
    C₂  = sim.bnd.weights[3]
    C₂ᵀ = sim.bnd.weights[4]

    𝕀1 = sparse(complex(1.,0)I, sim.dis.N[1], sim.dis.N[1])
    𝕀2 = sparse(complex(1.,0)I, sim.dis.N[2], sim.dis.N[2])

    return (
    (𝕀2 ⊗ ∇₁²) + (∇₂² ⊗ 𝕀1) +
    sparse(exp(+1im*ka*sim.lat.a)*I,size(C₁))*C₁ +
    sparse(exp(-1im*ka*sim.lat.a)*I,size(C₁ᵀ))*C₁ᵀ +
    sparse(exp(+1im*kb*sim.lat.b)*I,size(C₂))*C₂ +
    sparse(exp(-1im*kb*sim.lat.b)*I,size(C₂ᵀ))*C₂ᵀ
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

        C₁ = spzeros(Float64,prod(sim.dis.N),prod(sim.dis.N))
        C₁ᵀ = spzeros(Float64,prod(sim.dis.N),prod(sim.dis.N))
        C₂ = spzeros(Float64,prod(sim.dis.N),prod(sim.dis.N))
        C₂ᵀ = spzeros(Float64,prod(sim.dis.N),prod(sim.dis.N))

        if :p ∈ bc[:,1] && bc[1,1] == bc[2,1]
            I, J, V =  periodic_boundary_weights(sim, 1)
            C₁  = sparse(I, J, V/dx², prod(sim.dis.N), prod(sim.dis.N))
            C₁ᵀ = sparse(J, I, V/dx², prod(sim.dis.N), prod(sim.dis.N))
        elseif :p ∈ bc[:,1]
            throw(ArgumentError("only one boundary of dimension 1 is periodic, must be both or none"))
        end

        if :p ∈ bc[:,2] && bc[1,2] == bc[2,2]
            I, J, V = periodic_boundary_weights(sim, 2)
            C₂  = sparse(I, J, V/dy², prod(sim.dis.N), prod(sim.dis.N))
            C₂ᵀ = sparse(J, I, V/dy², prod(sim.dis.N), prod(sim.dis.N))
        elseif :p ∈ bc[:,2]
            throw(ArgumentError("only one boundary of dimension 2 is periodic, must be both or none"))
        end

        sim.bnd.weights[1] = C₁
        sim.bnd.weights[2] = C₁ᵀ
        sim.bnd.weights[3] = C₂
        sim.bnd.weights[4] = C₂ᵀ
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
    x = sim.dis.x[1]; y = sim.dis.x[2]
    dx = sim.dis.dx[1]; dy = sim.dis.dx[2]

    if dim == 1
        XY = bravais_coordinates_unit_cell.(x[1]-dx, y, Ref(lattice))
    else
        XY = bravais_coordinates_unit_cell.(x, y[1]-dy, Ref(lattice))
    end

    X = Array{Float64}(undef, N[mod1(dim+1,2)])
    Y = Array{Float64}(undef, N[mod1(dim+1,2)])
    for i ∈ eachindex(XY)
        X[i] = XY[i][1]
        Y[i] = XY[i][2]
    end

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

    I = Array{Int}(undef,4N[mod1(dim+1,2)])
    J = Array{Int}(undef,4N[mod1(dim+1,2)])
    V = Array{Float64}(undef,4N[mod1(dim+1,2)])

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

        j[:], k[:], l[:] = findnz( sparse(t, u, v, N[2], N[2]) ⊗ sparse(q, r, s, N[1], N[1]) )

        I[(4(i-1)+1):(4(i-1)+4)] = j
        J[(4(i-1)+1):(4(i-1)+4)] = k
        V[(4(i-1)+1):(4(i-1)+4)] = l

        # next!(pg)
    end
    return I, J, V
end


"""
    s₁, s₂ = σ(sim)

conductivity for absorbing layer (PML or not) in dimensions 1 and 2.
"""
function σ(sim::Simulation)

    α = Array{ComplexF64}(undef,2,2)
    for i ∈ CartesianIndices(α)
        if sim.bnd.bl[i] !== :none
            α[i] = -(1/4)*(float(POWER_LAW)+1)*exp(complex(0,SCALING_ANGLE::Float64))*log(EXTINCTION::Float64)/(sim.bnd.bl_depth[i]^(float(POWER_LAW)+1))
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
