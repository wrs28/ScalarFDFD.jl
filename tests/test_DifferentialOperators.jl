push!(LOAD_PATH,"/Users/wrs/wrs julia/ScalarFDFD/src")
using BenchmarkTools
using Bravais
using BoundaryConditions
using CoordinateSystems

using Revise
using DifferentialOperators

##

N = (1000,1001)
Δ1 = [0 0; 1 1]
Δ2 = ((0,0),(1,1))
Δ3 = [0 1 0 1]
Δ4 = (0,1,0,1)
# bcs = (DirichletBC{1,1}(),PeriodicBC{2,2}(BravaisLattice(a=1,b=1,β=π/3)),RobinBC{2,1}(1,sin),MatchedBC{1,2}([1],[2]))
bcs = ezbc(:d;lattice=BravaisLattice(a=1,b=1))
# BCS = ((bcs[1],bcs[2]),(bcs[3],bcs[4]))
OD = OperatorDefinition{2,1,Cartesian}(N,Δ4,bcs)

# @btime OperatorDefinition{2,1,Polar}(N,Δ1,bcs,ones(ComplexF64,N[1]+1,N[2]+1))
# @code_warntype OperatorDefinition{2,1,Polar}(N,Δ1,bcs,ones(ComplexF64,N[1]+1,N[2]+1))
#
# @btime OperatorDefinition{2,1,Polar}(N,Δ2,bcs)
# @code_warntype OperatorDefinition{2,1,Polar}(N,Δ4,bcs)

# @btime DifferentialOperators.RowColVal((10,11))
# @code_warntype DifferentialOperators.RowColVal((10,11))

OC = DifferentialOperators.OperatorConstructor(OD)
@btime DifferentialOperators.OperatorConstructor(OD)
@code_warntype DifferentialOperators.OperatorConstructor(OD)
