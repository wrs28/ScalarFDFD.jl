### DOMAIN STRUCT
Base.show(io::IO, dom::Domain) = begin
    print(io, typeof(dom), ": \n")
    print(io, "\tDomain type: ", dom.domain_type, "\n")
    if dom.which_asymptote !== :none
        print(io, "\t\twaveguide: ", dom.which_waveguide,"\n")
        print(io, "\t\tasymptotic region: ", dom.which_asymptote, "\n")
    end
    print(io, "\tdomain function: ", dom.is_in_domain, "\n")
    print(IOContext(io, :typeinfo => Array{Float64}), "\tdomain params: ", dom.domain_params, "\n")
    print(IOContext(io, :sub=>true), "\tdomain lattice: ", dom.lattice, "\n")
    print(IOContext(io, :typeinfo => Array{Float64}), "\tbackground Re{n}: ", dom.n₁, "\n")
    print(IOContext(io, :typeinfo => Array{Float64}), "\tbackground Im{n}: ", dom.n₂, "\n")
    print(IOContext(io, :typeinfo => Array{Float64}), "\tbackground F: ", dom.F, "\n")
    print(IOContext(io, :typeinfo => Array{Float64}), "\tRe{n} by region: ", dom.n₁_val[dom.n₁_idx], "\n")
    print(IOContext(io, :typeinfo => Array{Float64}), "\tIm{n} by region: ", dom.n₂_val[dom.n₂_idx], "\n")
    print(IOContext(io, :typeinfo => Array{Float64}), "\tF by region: ", dom.F_val[dom.F_idx], "\n")
    print(IOContext(io, :typeinfo => Array{Function}), "\tregion functions: ", dom.is_in_region, "\n")
    print(IOContext(io, :typeinfo => Array{Array{Float64,1},1}), "\tregion params: ", dom.region_params)
end


### SYSTEM STRUCT
Base.show(io::IO, sys::System) = begin
    if !get(io, :sub, false)
        print(io, typeof(sys), " with ", length(sys.domains), " domains: \n")
    end
    domain_string = [["\tdomain ", i, " type: ", sys.domains[i].domain_type, "\n"] for i ∈ eachindex(sys.domains)]
    for i ∈ eachindex(sys.domains)
        print(io, domain_string[i]...)
    end
    for i ∈ eachindex(sys.waveguides)
        if i==1
            print(io,"\n")
        end
        print(IOContext(io, :typeinfo => Array{Int}), "\twaveguide ", sys.waveguides[i], " domains: ", findall(sys.waveguides[i].==[sys.domains[j].which_waveguide for j ∈ eachindex(sys.domains)]), "\n")
        if i==length(sys.waveguides)
            print(io,"\n")
        end
    end
    # print(IOContext(io, :typeinfo => Array{Int}), "\tdomain: ", sys.domain_by_region, "\n")
    # print(IOContext(io, :typeinfo => Array{ComplexF64}), "\tn: ", sys.n_by_region, "\n")
    # print(IOContext(io, :typeinfo => Array{ComplexF64}), "\tε: ", sys.ε_by_region, "\n")
    # print(IOContext(io, :typeinfo => Array{Float64}), "\tF: ", sys.F_by_region)
end


### BOUNDARY STRUCT
Base.show(io::IO, bnd::Boundary) = begin
    if !get(io, :sub, false)
        print(io, typeof(bnd), ": \n")
    end

    print(io, "\t\t-------------------------\n")

    print(io, "\tbound.\t|\t   ", fmt("6s",bnd.bc[2,2]), "\t|\n")
    print(io, "\tcond-\t|   ",fmt("6s",bnd.bc[1,1]), "\t   ", fmt("5s",bnd.bc[2,1]), "|\n")
    print(io, "\t ition\t|\t   ", fmt("6s",bnd.bc[1,2]), "\t|\n")

    print(io, "\t\t-------------------------\n")

    print(io, "\t\t|\t", fmt("+5.3f",bnd.∂Ω[2,2]), "\t\t|\n")
    print(io, "\t∂Ω\t| ", fmt("+5.3f",bnd.∂Ω[1,1]), "\t", fmt("+5.3f",bnd.∂Ω[2,1]), "\t|\n")
    print(io, "\t\t|\t", fmt("+5.3f",bnd.∂Ω[1,2]), "\t\t|\n")

    print(io, "\t\t-------------------------\n")

    print(io, "\t\t|\t", fmt("+5.3f",bnd.∂Ω_tr[2,2]), "\t\t|\n")
    print(io, "\t∂Ω_tr\t| ",fmt("+5.3f",bnd.∂Ω_tr[1,1]), "\t", fmt("+5.3f",bnd.∂Ω_tr[2,1]), "\t|\n")
    print(io, "\t\t|\t", fmt("+5.3f",bnd.∂Ω_tr[1,2]), "\t\t|\n")

    print(io, "\t\t-------------------------\n")

    print(io, "\tbound.\t|\t", fmt("8s",bnd.bl[2,2]), "\t|\n")
    print(io, "\tlayer\t|",fmt("8s",bnd.bl[1,1]), "\t", fmt("8s",bnd.bl[2,1]), "|\n")
    print(io, "\t\t|\t", fmt("8s",bnd.bl[1,2]), "\t|\n")

    print(io, "\t\t-------------------------\n")

    print(io, "\tbound.\t|\t", fmt("+5.3f",bnd.bl_depth[2,2]), "\t\t|\n")
    print(io, "\tlayer\t| ",fmt("+5.3f",bnd.bl_depth[1,1]), "\t", fmt("+5.3f",bnd.bl_depth[2,1]), "\t|\n")
    print(io, "\tdepth\t|\t", fmt("+5.3f",bnd.bl_depth[1,2]), "\t\t|\n")

    print(io, "\t\t-------------------------\n")
end


### DISCRETIZATION STRUCT
Base.show(io::IO, dis::Discretization) = begin
    if !get(io, :sub, false)
        print(io, typeof(dis), ": \n")
    end
    print(io, "\tN: ", dis.N, "\n",
    # "\ttruncated N: ", dis.N_tr, "\n",
    "\tsub-pixel number: ", dis.sub_pixel_num, "\n",
    "\tdx: ", dis.dx)
end


### CHANNEL STRUCT
Base.show(io::IO, chn::Channels) = begin
    if get(io, :indented, false)
        print(io, "\t\twaveguide: ", chn.waveguide, "\n",
        "\t\tquantum number: ", chn.quantum_number)
    else
        print(io, typeof(chn), ": \n")
        print(io, "\twaveguide: ", chn.waveguide, "\n",
        "\tquantum number: ", chn.quantum_number, "\n")
        if !isempty(chn.dispersion)
            dispersion_plot = lineplot(x->chn.dispersion[1](x), 0:.01:2π, canvas=BlockCanvas, color=:yellow, title="Dispersion", name="channel", ylim=[chn.gaps[1][1]-.1, chn.gaps[1][2]+.1]);
            lineplot!(dispersion_plot, chn.gaps[1][1], 0, color=:white);
            lineplot!(dispersion_plot, chn.gaps[1][2], 0, color=:white);
            println(dispersion_plot)
        end
    end
end


### SCATTERING STRUCT
Base.show(io::IO, sct::Scattering) = begin
    plot_flag=false
    if !get(io, :sub, false)
        plot_flag = true
        print(io, typeof(sct), " with ", length(sct.channels), " channels:\n")
    end
    temp = [["\n", sct.channels[i]] for i ∈ eachindex(sct.channels)]
    if plot_flag && !isempty(sct.channels) && !isempty(sct.channels[1].dispersion)
        dispersion_plot = lineplot(x->sct.channels[1].dispersion[1](x), 0:.01:2π, canvas=BlockCanvas, title="Dispersion", name="channel 1", ylim=[sct.channels[1].gaps[1][1]-.1, sct.channels[1].gaps[1][2]+.1])
    end
    for i ∈ eachindex(temp)
        if i>1
            print(io, "\n")
        end
        print(io, "\tChannel ", i, ":")
        print(IOContext(io::IO, :indented => true), temp[i]...)
        if plot_flag
            if !isempty(sct.channels[i].dispersion)
                if i>1
                    lineplot!(dispersion_plot,x->sct.channels[i].dispersion[1](x), 0:.01:2π, name="channel $i");
                end
                lineplot!(dispersion_plot, sct.channels[i].gaps[1][1], 0, color=:white);
                lineplot!(dispersion_plot, sct.channels[i].gaps[1][2], 0, color=:white);
            end
        end
    end
    try
        dispersion_plot
        print(io,"\n\n")
        println(dispersion_plot)
    catch
        nothing
    end
end


### TWO-LEVEL-SYSTEM STRUCT
Base.show(io::IO, tls::TwoLevelSystem) = begin
    if !get(io, :sub, false)
        print(io, typeof(tls), ":\n")
    end
    print(io, "\tD₀: ", tls.D₀, "\n",
    "\tω₀: ", tls.ω₀, "\n",
    "\tγ⟂: ", tls.γp)
end


### SIMULATION STRUCT
Base.show(io::IO, sim::Simulation) = begin
    print(IOContext(io, :sub=>true),
    typeof(sim), ": \n\n",
    "sys: ", sim.sys, "\n\n",
    "bnd: ", sim.bnd, "\n\n",
    "dis: ", sim.dis, "\n\n",
    "sct: ", sim.sct, "\n\n",
    "tls: ", sim.tls)
    if [sim.lat.a, sim.lat.b] !== [Inf,Inf]
        print(io, "\n\nlattice: ", sim.lat)
    end
end

### BRAVAIS LATTICE STRUCT
Base.show(io::IO, bvl::Bravais) = begin
    if !get(io, :sub, false)
        print(io, typeof(bvl), ": \n",
        "\tprimitive vector 1: ", fmt("3.2f",bvl.a), ", ∠", fmt("3.2f",(mod2pi(bvl.α))*180/pi), "°\n",
        "\tprimitive vector 2: ", fmt("3.2f",bvl.b), ", ∠", fmt("3.2f",(mod2pi(bvl.β))*180/pi), "°\n",
        "\torigin: (", fmt("3.2f",bvl.x0), ", ", fmt("3.2f",bvl.y0), ")")
    else
        print(io,
        "\n\t\tprimitive vector 1: ", fmt("3.2f",bvl.a), ", ∠", fmt("3.2f",(mod2pi(bvl.α))*180/pi), "°\n",
        "\t\tprimitive vector 2: ", fmt("3.2f",bvl.b), ", ∠", fmt("3.2f",(mod2pi(bvl.β))*180/pi), "°\n",
        "\t\torigin: (", fmt("3.2f",bvl.x0), ", ", fmt("3.2f",bvl.y0), ")")
    end
end