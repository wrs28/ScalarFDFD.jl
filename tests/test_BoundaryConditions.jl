push!(LOAD_PATH,"/Users/wrs/wrs julia/ScalarFDFD/src")
using BenchmarkTools
using Bravais
using Revise
using BoundaryConditions

@btime PML{1,2}(1)
@code_warntype PML{1,2}(1)

@btime cPML{2,1}(1)
@code_warntype cPML{2,1}(1)

@btime noBL{2,2}()
@code_warntype noBL{2,2}()

@btime DirichletBC{2,1}()
@code_warntype DirichletBC{2,1}()

@btime DirichletBC{2,1}(0)
@code_warntype DirichletBC{2,1}(0)

@btime DirichletBC{2,1}(sin)
@code_warntype DirichletBC{2,1}(sin)

@btime DirichletBC{2,1}(1,(10,11))
@code_warntype DirichletBC{2,1}(1,(10,11))

@btime DirichletBC{2,1}(sin,(10,11);dx=(.1,.1))
@code_warntype DirichletBC{2,1}(sin,(10,11);dx=(.1,.1))

@btime DirichletBC{1,2}(sin.(1:100),(100,101);dx=(.1,.1))
@code_warntype DirichletBC{1,2}(sin.(1:100),(100,101);dx=(.1,.1))

dbc = DirichletBC{1,1}()
@btime dbc((100,101))
@code_warntype dbc((100,101))

dbc = DirichletBC{1,1}(sin)
@btime dbc((100,101);dx=(.1,.1))
@code_warntype dbc((100,101);dx=(.1,.1))

DBC = dbc((100,101);dx=(.1,.1))
@btime DBC((10,11);dx=(.1,.1))
@code_warntype DBC((10,11);dx=(.1,.1))


@btime NeumannBC{2,1}()
@code_warntype NeumannBC{2,1}()

@btime NeumannBC{2,1}(0)
@code_warntype NeumannBC{2,1}(0)

@btime NeumannBC{2,1}(sin)
@code_warntype NeumannBC{2,1}(sin)

@btime NeumannBC{2,1}((10,11))
@code_warntype NeumannBC{2,1}((10,11))

@btime NeumannBC{2,1}(sin,(10,11);dx=(.1,.1))
@code_warntype NeumannBC{2,1}(sin,(10,11);dx=(.1,.1))

@btime NeumannBC{1,2}(sin.(1:100),(100,101);dx=(.1,.1))
@code_warntype NeumannBC{1,2}(sin.(1:100),(100,101);dx=(.1,.1))

nbc = NeumannBC{1,1}()
@btime nbc((100,101);dx=(.1,.2))
@code_warntype nbc((100,101);dx=(.1,.2))

nbc = NeumannBC{1,1}(sin)
@btime nbc((100,101);dx=(.1,.1))
@code_warntype nbc((100,101);dx=(.1,.1))

NBC = nbc((100,101);dx=(.1,.1))
@btime NBC((10,11);dx=(.1,.1))
@code_warntype NBC((10,11);dx=(.1,.1))


@btime RobinBC{2,1}(1,.2)
@code_warntype RobinBC{2,1}(1,.2)

@btime RobinBC{2,1}(1,sin)
@code_warntype RobinBC{2,1}(1,sin)

@btime RobinBC{2,1}(sin,1)
@code_warntype RobinBC{2,1}(1,sin)

@btime RobinBC{2,1}(sin,cos)
@code_warntype RobinBC{2,1}(cos,sin)

@btime RobinBC{2,1}(1,1,1,(10,11);dx=(.1,.1))
@code_warntype RobinBC{2,1}(1,1,1,(10,11);dx=(.1,.1))

@btime RobinBC{2,1}(sin,1,1,(10,11);dx=(.1,.1))
@code_warntype RobinBC{2,1}(sin,1,1,(10,11);dx=(.1,.1))

@btime RobinBC{1,2}(1,sin,1,(10,11);dx=(.1,.1))
@code_warntype RobinBC{2,2}(1,sin,1,(10,11);dx=(.1,.1))

@btime RobinBC{2,1}(1,1,cos,(10,11);dx=(.1,.1))
@code_warntype RobinBC{2,1}(1,1,cos,(10,11);dx=(.1,.1))

@btime RobinBC{2,1}(sin,cos,1,(10,11);dx=(.1,.1))
@code_warntype RobinBC{2,1}(sin,cos,1,(10,11);dx=(.1,.1))

@btime RobinBC{2,1}(sin,1,cos,(10,11);dx=(.1,.1))

@btime RobinBC{2,1}(1,exp,cos,(10,11);dx=(.1,.1))
@code_warntype RobinBC{2,1}(1,exp,cos,(10,11);dx=(.1,.1))

@btime RobinBC{2,2}(sin,exp,cos,(10,11);dx=(.1,.1))
@code_warntype RobinBC{2,2}(sin,exp,cos,(10,11);dx=(.1,.1))

@btime RobinBC{2,2}(rand(10,11),rand(10,11),rand(10,11),(10,11);dx=(.1,.1))
@code_warntype RobinBC{2,2}(rand(10,11),rand(10,11),rand(10,11),(10,11);dx=(.1,.1))

rbc = RobinBC{2,1}(sin,cos,exp)
RBC = rbc((10,11);dx=(.1,.1))
@btime rbc((10,11);dx=(.1,.1))
@code_warntype rbc((10,11);dx=(.1,.1))


lat = BravaisLattice(a=1,b=1,α=.2,β=π/3)

@btime PeriodicBC{2,2}(lat)
@code_warntype PeriodicBC{2,2}(lat)

@btime PeriodicBC{1,1}(lat,(10,11))
@code_warntype PeriodicBC{1,1}(lat,(10,11))

pbc = PeriodicBC{2,1}(lat)
PBC = pbc((10,11))
@btime pbc((10,11))
@code_warntype pbc((10,11))

@btime PBC((10,11))
@code_warntype PBC((10,11))

@btime BoundaryConditions.get_bc_type(dbc)
@code_warntype BoundaryConditions.get_bc_type(dbc)
@btime BoundaryConditions.get_bc_type(nbc)
@code_warntype BoundaryConditions.get_bc_type(nbc)
@btime BoundaryConditions.get_bc_type(rbc)
@code_warntype BoundaryConditions.get_bc_type(rbc)

bc = (DirichletBC{1,1}(),NeumannBC{2,2}(),RobinBC{2,1}(1,sin),MatchedBC{1,2}([1],[2]))
@btime BoundaryConditions.get_dim_side(bc)
@code_warntype BoundaryConditions.get_dim_side(bc)

@btime BoundaryConditions.get_dim(bc)
@code_warntype BoundaryConditions.get_dim(bc)

@btime BoundaryConditions.get_side(bc)
@code_warntype BoundaryConditions.get_side(bc)

@btime BoundaryConditions.reorder_side(bc)
@code_warntype BoundaryConditions.reorder_side(bc)

@btime BoundaryConditions.reorder_dim(bc)
@code_warntype BoundaryConditions.reorder_dim(bc)

@btime BoundaryConditions.reorder(bc)
@code_warntype BoundaryConditions.reorder(bc)

@btime BoundaryConditions.apply_args_to_bc(bc;N=(10,11),dx=(.1,.1),coordinate_system=Polar())
@code_warntype BoundaryConditions.apply_args_to_bc(bc;N=(10,11),dx=(.1,.1),coordinate_system=Polar())

@btime ezbc(:d,:n)
@btime ezbc(:d,[:n,:d])
@btime ezbc([:o,:d],:p;lattice=lat,outgoing_qns=[1],incoming_qns=[2])
@btime ezbc([:o,:d],[:p,:r];lattice=lat,outgoing_qns=[1],incoming_qns=[2],α=1,β=2)

@btime MatchedBC{1,2}()
@code_warntype MatchedBC{1,2}()

@btime MatchedBC{1,2,Polar}()
@code_warntype MatchedBC{1,2,Polar}()

@btime MatchedBC{1,2,Polar}()
@code_warntype MatchedBC{1,2,Polar}()

@btime MatchedBC{1,2,Polar}([1],[2])
@code_warntype MatchedBC{1,2,Polar}([1],[2])

@btime MatchedBC{1,2,Polar}([1],[2],(10,11);dx=(.1,.1))
@code_warntype MatchedBC{1,2,Polar}([1],[2],(10,11);dx=(.1,.1))

mbc = MatchedBC{1,2,Polar}([1],[2])
MBC = mbc((10,11);dx=(.1,.1),coordinate_system=Polar())
@btime mbc((10,11);dx=(.1,.1),coordinate_system=Polar())
@code_warntype mbc((10,11);dx=(.1,.1),coordinate_system=Polar())
@btime MBC((10,11);dx=(.1,.1),coordinate_system=Polar())
@code_warntype MBC((10,11);dx=(.1,.1),coordinate_system=Polar())
