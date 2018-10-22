################################################################################
# REGIONS STRUCTURE
################################################################################
"""
    regions = Regions()
    regions = Regions(functions, params, n1, n2, F)
    regions = Regions(domain)

`Regions()` initializes empty region

`Regions(functions, params, n1, n2, F)` creates region with specified parameters

`Regions(domain)` makes a region out of background default domain parameters (i.e. n1, n2, etc)
"""
struct Regions
    functions::Array{Function,1}
    params::Array{Array{Float64,1},1}
    n₁::Array{Float64,1}
    n₂::Array{Float64,1}
    F::Array{Float64,1}

    function Regions()
        new(Function[], Array{Float64}[], Float64[], Float64[], Float64[])
    end
    function Regions(functions, params, n₁, n₂, F)
        new(functions, float.(params), float.(n₁), float.(n₂), float.(F))
    end
    function Regions(domain::Domain)
        return Regions(domain.is_in_region, domain.region_params, domain.n₁_val, domain.n₂_val, domain.F_val)
    end
end
Base.show(io::IO, rgn::Regions) = begin
    print(io,
    typeof(rgn), " with ", length(rgn.functions), " regions:")
    temp1 = [["\n\t\tfunction: ", rgn.functions[i]] for i ∈ eachindex(rgn.functions)]
    temp2 = [["\n\t\tparams: ", rgn.params[i]] for i ∈ eachindex(rgn.params)]
    temp3 = [["\n\t\tindex n: ", complex(rgn.n₁[i],rgn.n₂[i])] for i ∈ eachindex(rgn.n₁)]
    temp4 = [["\n\t\tpump F: ", rgn.F[i]] for i ∈ eachindex(rgn.F)]
    for i ∈ eachindex(temp1)
        print(io, "\n\tregion ", i)
        print(io, temp1[i]...)
        print(io, temp2[i]...)
        print(io, temp3[i]...)
        print(io, temp4[i]...)
    end
end


################################################################################
# DOMAIN BUILDING
################################################################################
"""
    domain = build_domain(regions; is_in_domain=whole_domain, domain_params=[],
                n₁=1, n₂=0, F=0, domain_type=:generic, lattice=Bravais())

Builds a domain out of `regions` (::Regions), with background parameters
`is_in_domain`, `domain_params`, etc.
"""
function build_domain(regions::Regions=Regions(); is_in_domain=whole_domain, domain_params=[],
                        n₁=1, n₂=0, F=0, domain_type=:generic, lattice=Bravais(),
                        which_asymptote=:none, which_waveguide=0)
    return Domain(; is_in_domain=is_in_domain, n₁=n₁, n₂=n₂, F=F,
                    domain_params=domain_params, domain_type=domain_type,
                    is_in_region=regions.functions, region_params=regions.params,
                    n₁_val=regions.n₁, n₁_idx=1:length(regions.n₁),
                    n₂_val=regions.n₂, n₂_idx=1:length(regions.n₂),
                    F_val=regions.F, F_idx=1:length(regions.n₂),
                    lattice=lattice, which_asymptote=which_asymptote, which_waveguide=which_waveguide )
end


################################################################################
# GENERIC REGION CONSTRUCTION
################################################################################
"""
    new_domain = add_regions(new_function, old_domain)

    new_regions = add_regions(regions, old_regions)

add a new region either to domain or to a list of regions
"""
function add_regions(new_regions::Regions, domain::Domain)
    return build_domain(add_regions(new_regions, Regions(domain);
            is_in_domain=domain.is_in_domain, domain_params=domain.domain_params,
            n₁=domain.n₁, n₂=domain.n₂, F=domain.F, domain_type=domain.domain_type,
            lattice=domain.lattice))
end
function add_regions(new_regions::Regions, old_regions::Regions=Regions())
    return Regions(
                new_region.functions, new_regions.params,
                new_region.n₁, new_regions.n₂, new_regions.F,
                old_regions
                )
end


"""
    new_domain = add_region(new_function, new_params, old_domain, new_n₁=1, new_n₂=0, new_F=0)

    new_regions = add_region(new_function, new_params, new_n₁=1, new_n₂=0, new_F=0, old_regions=Regions())

add a new region either to `old_domain` or to `old_regions`.

`old_regions` will initialize to an empty region if not provided.
"""
function add_region(new_region_function::Function, new_region_params, domain::Domain,
                    new_region_n₁=1, new_region_n₂=0, new_region_F=0)
    old_regions = Regions(domain)
    regions = add_region(new_region_function, new_region_params, new_region_n₁, new_region_n₂, new_region_F, old_regions)
    return build_domain(regions; is_in_domain=domain.is_in_domain, domain_params=domain.domain_params,
                            n₁=domain.n₁, n₂=domain.n₂, F=domain.F, domain_type=domain.domain_type, lattice=domain.lattice)
end

function add_region(new_region_function::Function, new_region_params,
                    new_region_n₁=1, new_region_n₂=0, new_region_F=0,
                    old_regions::Regions=Regions() )
    return Regions(
                vcat(new_region_function, old_regions.functions),
                vcat([new_region_params], old_regions.params),
                vcat(new_region_n₁, old_regions.n₁),
                vcat(new_region_n₂, old_regions.n₂),
                vcat(new_region_F, old_regions.F)
                )
end
