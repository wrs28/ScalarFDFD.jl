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


"""
    isDirichlet(bc)
"""
function isDirichlet(bc::Symbol)
    validNames = [:d, :D, :dirichlet, :Dirichlet, :hard, :h]
    return bc ∈ validNames
end


"""
    isNeumann(bc)
"""
function isNeumann(bc::Symbol)
    validNames = [:d, :D, :dirichlet, :Dirichlet, :hard, :h]
    return bc ∈ validNames
end


"""
    isOpen(bc)
"""
function isOpen(bc::Symbol)
    validNames = [:o, :open, :Open]
    return bc ∈ validNames
end


"""
    isPeriodic(bc)
"""
function isPeriodic(bc::Symbol)
    validNames = [:p, :periodic, :Periodic, :bloch, :Bloch]
    return bc ∈ validNames
end


"""
    isPMLout(bl)
"""
function isPMLout(bl::Symbol)
    validNames = [:pml_out, :PML_OUT, :PML_out, :pml, :PML]
    return bl ∈ validNames
end


"""
    isPMLin(bl)
"""
function isPMLin(bl::Symbol)
    validNames = [:pml_in, :PML_IN, :PML_in, :pml_conj, :PML_conj]
    return bl ∈ validNames
end


"""
    isABSout(bl)
"""
function isABSout(bl::Symbol)
    validNames = [:abs_out, :ABS_OUT, :ABS_out, :abs, :ABS, :amp, :AMP]
    return bl ∈ validNames
end


"""
    isABSin(bl)
"""
function isABSin(bl::Symbol)
    validNames = [:abs_in, :ABS_IN, :ABS_in, :abs_conj, :ABS_conj]
    return bl ∈ validNames
end


"""
    isNone(bl)
"""
function isNone(bl::Symbol)
    validNames = [:none, :nothing, :empty]
    return bl ∈ validNames
end
