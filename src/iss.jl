"""
    iswaveguide(domain::Domain)
"""
function iswaveguide(domain::Domain)
    return domain.domain_type ∈ [:planar_waveguide, :planar_waveguide_background, :halfspace, :pc_waveguide, :pc_waveguide_background]
end


"""
    isbulkwaveguide(domain::Domain))
"""
function isbulkwaveguide(domain::Domain)
    return domain.domain_type ∈ [:bulk_planar_waveguide_x, :bulk_planar_waveguide_y, :bulk_pc_waveguide_x, :bulk_pc_waveguide_y]
end


"""
    isbackground(domain::Domain)
"""
function isbackground(domain)
    return domain.domain_type ∈ [:background, :planar_waveguide_background, :pc_waveguide_background]
end


"""
    isdefect(domain::Domain))
"""
function isdefect(domain::Domain)
    return domain.domain_type ∈ [:defect, :site_defect, :line_defect_simple, :line_defect_periodic, :pc_waveguide]
end


"""
    ispc(dom::Domain)
"""
function ispc(dom::Domain)
    return dom.domain_type ∈ [:pc, :pc_waveguide, :pc_waveguide_background]
end


"""
    isplanar(dom::Domain)
"""
function isplanar(dom::Domain)
    return dom.domain_type ∈ [:planar_waveguide, :planar_waveguide_background, :bulk_planar_waveguide_x, :bulk_planar_waveguide_y]
end


"""
    ishalfspace(dom::Domain)
"""
function ishalfspace(dom::Domain)
    return dom.domain_type ∈ [:halfspace_waveguide]
end
