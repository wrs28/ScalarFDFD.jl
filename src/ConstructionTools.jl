module ConstructionTools

using Shapes,
SimulationDefinition

export Subdomain

"""
    subdomains = Subdomains()
    subdomains = Subdomains(functions, params, type::Shape, ε::Function, F::Function)
    subdomains = Subdomains(domain)

`Regions()` initializes empty region

`Regions(functions, params, n1, n2, F)` creates region with specified parameters

`Regions(domain)` makes a region out of background default domain parameters (i.e. n1, n2, etc)
"""
struct Subdomains{Tp}
    is_in_subdomain::Array{Function,1}
    params::Array{Tp,1}
    type::Array{Symbol,1}
    ε::Array{Function,1}
    F::Array{Function,1}

    function Subdomains()
        new{Dict{Symbol,Float64},Nothing}(Function[], Dict{Symbol,Float64}[], Symbol[], Function[], Function[])
    end
    function Subdomains(is_in_subdomain, params, type, ε, F)
        new(is_in_subdomain, params, type, ε, F)
    end
    function Subdomains(dom::Tdom) where Tdom<:Domain
        return Subdomains(dom.is_in_subdomain, dom.subdomain_params, dom.subdomain_type, dom.subdomain_ε, dom.subdomain_F)
    end
end
Base.show(io::IO, sbd::Subdomains) = begin
    print(io,
    typeof(sbd), " with ", length(sbd.is_in_subdomain), " regions:")
    temp1 = [["\n\t\tfunction: ", sbd.params[i]] for i ∈ eachindex(sbd.params)]
    temp2 = [["\n\t\tparams: ", sbd.type[i]] for i ∈ eachindex(sbd.type)]
    for i ∈ eachindex(temp1)
        print(io, "\n\tregion ", i)
        print(io, temp1[i]...)
        print(io, temp2[i]...)
    end
end


"""
    domain = build_domain(subdomains; is_in_domain=whole_domain, domain_params=[],
                n₁=1, n₂=0, F=0, domain_type=:generic, lattice=Bravais())

Builds a domain out of `regions` (::Regions), with background parameters
`is_in_domain`, `domain_params`, etc.
"""
function build_domain(subomains::Tsd=Subdomains();
    is_in_domain::Function=whole_domain,
    domain_params::Td=Dict(:n₁ => 1, :n₂ => 0, :F => 0),
    domain_type::Symbol=:background,
    domain_ε::Function = piecewise_constant_ε,
    domain_F::Function = piecewise_constant_F,
    lattice::BravaisLattice = BravaisLattice(),
    which_asymptote::Symbol = :none,
    which_waveguide::Int = 0
    ) where Tsd<:Subdomain where Td<:Dict

    return Domain(; is_in_domain=is_in_domain,
                    domain_params=domain_params,
                    domain_type=domain_type,
                    domain_ε=domain_ε,
                    domain_F=domain_F,
                    is_in_subdomain = subdomain.is_in_subdomain,
                    subdomain_params = subdomain.params,
                    subdomain_type = subdomain.type,
                    subdomain_ε = subdomain.ε,
                    subdomain_F = subdomain.F,
                    lattice=lattice,
                    which_asymptote=which_asymptote,
                    which_waveguide=which_waveguide
                    )
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

end #module
