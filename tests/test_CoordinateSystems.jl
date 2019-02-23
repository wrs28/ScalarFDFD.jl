push!(LOAD_PATH,"/Users/wrs/wrs julia/ScalarFDFD/src")
using BenchmarkTools
using Revise
using CoordinateSystems

cart = Cartesian()
pol = Polar()

@btime isCartesian(cart)
@code_warntype isCartesian(cart)

@btime isCartesian(pol)
@code_warntype isCartesian(pol)

@btime isPolar(pol)
@code_warntype isPolar(cart)

@btime isPolar(pol)
@code_warntype isPolar(pol)
