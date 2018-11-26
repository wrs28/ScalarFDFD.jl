################################################################################
#### RESONANT EIGENVALUE SOLVERS
################################################################################

"""
Linear Eigenvalue Solvers:
    freqs, ψ = eig_k(sim, k, nk=1; ka=0, kb=0, F=[1], ψ_init=[], is_linear=true)
    freqs, ψ = eig_k(sim, k, k_type::Symbol, nk=1; direction=[1,1], ka=0, kb=0, F=[1], ψ_init=[], is_linear=true)

Nonlinear Eigenvalue Solvers:
    freqs, ψ = eig_k(sim, k;  ka=0, kb=0, η_init=0, F=[1], ψ_init=[], is_linear=false)
    freqs, ψ = eig_k(sim, k, k_type::Symbol, nk=1, radii=(.1,.1); direction=[1,1], ka=0, kb=0, F=[1], ψ_init=[], Nq=100, is_linear=false)

Outputs:

- `freqs` = `freqs`[# of poles/zeros/etc]

- `ψ` = `ψ`[cavity size, # of poles]


Inputs (both linear and nonlinear):

- `k` is the frequency anchor, and those solutions closest to it are found

- `k_type` ∈ [:Z,:Zero,:z,:zero,:P,:Pole,:p,:pole,:UZR,:uzr,:U,:u]

- `nk` is the number of modes to be computed (default to 1), and may be ommitted.
When doing a nonlinear solve, the CF root-finding method does not take this argument.

- `ka`, `kb` are the bloch wavenumbers along each lattice vector
(doesn't matter unless at least one direction has periodic boundary conditions)

- `F` is an array which multiplies the bare pump profile defined in `sim`. This
is an easy way to include hole-burning.

- `ψ_init` is a vector which can seed the Arnoldi algorithm, in principle it is
useful for speeding up calculations if you know something near the solution. In
practice this has very little effect.


Nonlinear Inputs:

- `radii` are the axes of the elliptical contour

- `Nq` is the number of quadrature points along contour

- `η_init` is the initial guess for CF root-finding


Note: Set `D` or `F` to zero to get passive cavity.

Linear methods don't account for line-pulling, and only account for open boundaries
if either absorbing layers, PML's, or explicitly open boudnary conditions are used.
"""
function eig_k(sim::Simulation, k::Number; ka=0, kb=0, η_init=0,
            is_linear=true, k_avoid=[0], disp_opt=false, tol=.5, max_count=15, max_iter=50, args...)
    nk=1
    if is_linear
        k, ψ = eig_kl(sim, k, nk, ka, kb)
    else
        k, ψ, _ = eig_knl(sim, k, η_init, ka, kb, k_avoid, disp_opt, tol, max_count, max_iter)
    end
    return k::Array{ComplexF64,1}, ψ::Array{ComplexF64,2}
end # straight wrapper linear/CF root-finding


# straight wrapper linear/contour
function eig_k(sim::Simulation, k::Number, nk::Int, radii=(.1,.1); ka=0, kb=0, is_linear=true,
            r_min=.01, Nq=100, rank_tol=1e-8, parallel=nprocs()>1, args...)
    if is_linear
        k,ψ = eig_kl(sim, k, nk, ka, kb)
    else
        if parallel
            k = eig_knl(sim, k, nk, radii, ka, kb, Nq, rank_tol, r_min)
        else
            k = eig_knlp(sim, k, nk, radii, ka, kb, Nq, rank_tol, r_min)
        end
        ψ=fill(complex(NaN,NaN), prod(sim.dis.N), length(k))
    end
    return k::Array{ComplexF64,1}, ψ::Array{ComplexF64,2}
end


# k_type linear/CF root-finding
function eig_k(sim::Simulation, k::Number, k_type::Symbol; is_linear=true,
                direction::Array{Int,1}=[1,1], ka=0, kb=0, η_init=0, k_avoid=[0], disp_opt=false,
                tol=.5, max_count=15, max_iter=50, args...)

    bl_original = set_bl!(sim, k_type, direction)
    try
        k, ψ = eig_k(sim, k; is_linear=is_linear, ka=ka, kb=kb,
                        η_init=η_init, k_avoid=k_avoid, disp_opt=disp_opt, tol=tol,
                        max_count=max_count, max_iter=max_iter)
        return k, ψ
    finally
        reset_bl!(sim, bl_original)
    end
end


# k_type linear/contour
function eig_k(sim::Simulation, k::Number, k_type::Symbol, nk::Int, radii=(.1,.1);
            direction::Array{Int,1}=[1,1], is_linear=true, Nq=100,
            ka=0, kb=0, r_min=.01, rank_tol=1e-8, parallel=nprocs()>1, args...)

    bl_original = set_bl!(sim, k_type, direction)
    try
        k, ψ = eig_k(sim, k, nk, radii; is_linear=is_linear,
                        Nq=Nq, ka=ka, kb=kb, r_min=r_min, rank_tol=rank_tol, parallel=parallel)
        return k, ψ
    finally
        reset_bl!(sim, bl_original)
    end
end

################################################################################
#### CF EIGENVALUE SOLVER
################################################################################
"""
    u, η = eig_cf(sim, k, ncf=1; F=[1], η_init=0, u_init=[], ka=0, kb=0)
    u, η = eig_cf(sim, k, k_type::Symbol, ncf=1; F=[1], η_init=0, u_init=[], direction=[1,1], ka=0, kb=0)

CF eigenvalue solver.

- `k` is the frequency

- `k_type` ∈ [:Z,:Zero,:z,:zero,:P,:Pole,:p,:pole,:UZR,:uzr,:U,:u]

- `ncf` is the number of modes to be computed (default to 1), and may be ommitted.

- `ka`, `kb` are the bloch wavenumbers along each lattice vector
(doesn't matter unless at least one direction has periodic boundary conditions)

- `F` is an array which multiplies the bare pump profile defined in `sim`. This
is an easy way to include hole-burning.

- `u_init` is a vector which can seed the Arnoldi algorithm, in principle it is
useful for speeding up calculations if you know something near the solution. In
practice this has very little effect.

- `η_init` is the anchor CF eigenvalue, solutions are found which are closest to it.
"""
function eig_cf(sim::Simulation, k::Number, ncf::Int; η_init=0, ka=0, kb=0)

    η, u = eig_cf(sim, k, ncf, η_init, ka, kb)
    return η, u
end


function eig_cf(sim::Simulation, k, k_type::Symbol, ncf::Int=1; direction::Array{Int,1}=[1,1],
                    η_init=0, ka=0, kb=0)

    bl_original = set_bl!(sim, k_type, direction)
    try
        η, u = eig_cf(sim, k, ncf, η_init, ka, kb)
        return η, u
    finally
        reset_bl!(sim, bl_original)
    end
end


################################################################################
#### BOUNDARY FIXING FOR WRAPPERS
################################################################################
"""
    set_bl!(sim, k_type, direction)
"""
function set_bl!(sim::Simulation, k_type, direction=[])

    bl_original = deepcopy(sim.bnd.bl)

    if k_type ∈ [:Pole,:pole,:P,:p]
        for i ∈ eachindex(sim.bnd.bl)
            if sim.bnd.bl[i] == :pml_in
                sim.bnd.bl[i] = :pml_out
            elseif sim.bnd.bl[i] == :abs_in
                sim.bnd.bl[i] = :abs_out
            end
        end
    elseif k_type ∈ [:Zero,:zero,:Z,:z]
        for i ∈ eachindex(sim.bnd.bl)
            if sim.bnd.bl[i] == :pml_out
                sim.bnd.bl[i] = :pml_in
            elseif sim.bnd.bl[i] == :abs_out
                sim.bnd.bl[i] = :abs_in
            end
        end
    elseif k_type ∈ [:UZR,:uzr,:U,:u]
        for j ∈ 1:2
            if direction[j] == +1
                if sim.bnd.bl[1,j] == :pml_out
                    sim.bnd.bl[1,j] = :pml_in
                elseif sim.bnd.bl[1,j] == :abs_out
                    sim.bnd.bl[1,j] = :abs_in
                end
                if sim.bnd.bl[2,j] == :pml_in
                    sim.bnd.bl[2,j] = :pml_out
                elseif sim.bnd.bl[2,j] == :abs_in
                    sim.bnd.bl[2,j] = :abs_out
                end
            elseif direction[j] == -1
                if sim.bnd.bl[1,j] == :pml_in
                    sim.bnd.bl[1,j] = :pml_out
                elseif sim.bnd.bl[1,j] == :abs_in
                    sim.bnd.bl[1,j] = :abs_out
                end
                if sim.bnd.bl[2,j] == :pml_out
                    sim.bnd.bl[2,j] = :pml_in
                elseif sim.bnd.bl[2,j] == :abs_out
                    sim.bnd.bl[2,j] = :abs_in
                end
            else
                throw(ArgumentError("invalid direction $direction[j], should be ±1"))
            end
        end
    end
    return bl_original
end


"""
    reset_bl!(sim, k_type, direction)
"""
function reset_bl!(sim::Simulation, bl_original)
    @. sim.bnd.bl = bl_original
    return nothing
end
