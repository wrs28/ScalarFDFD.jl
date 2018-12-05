"""
    isWaveguide(domain::Domain)
"""
function isWaveguide(domain::Domain)
    return domain.domain_type ∈ [:planar_waveguide, :planar_waveguide_background, :halfspace, :pc_waveguide, :pc_waveguide_background]
end


"""
    isBulkWaveguide(domain::Domain))
"""
function isBulkWaveguide(domain::Domain)
    return domain.domain_type ∈ [:bulk_planar_waveguide_x, :bulk_planar_waveguide_y, :bulk_pc_waveguide_x, :bulk_pc_waveguide_y]
end


"""
    isBackground(domain::Domain)
"""
function isBackground(domain)
    return domain.domain_type ∈ [:background, :planar_waveguide_background, :pc_waveguide_background]
end


"""
    isDefect(domain::Domain))
"""
function isDefect(domain::Domain)
    return domain.domain_type ∈ [:defect, :site_defect, :line_defect_simple, :line_defect_periodic, :pc_waveguide]
end


"""
    isPC(dom::Domain)
"""
function isPC(dom::Domain)
    return dom.domain_type ∈ [:pc, :pc_waveguide, :pc_waveguide_background]
end


"""
    isPlanar(dom::Domain)
"""
function isPlanar(dom::Domain)
    return dom.domain_type ∈ [:planar_waveguide, :planar_waveguide_background, :bulk_planar_waveguide_x, :bulk_planar_waveguide_y]
end


"""
    isHalfSpace(dom::Domain)
"""
function isHalfSpace(dom::Domain)
    return dom.domain_type ∈ [:halfspace_waveguide]
end
