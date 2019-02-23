using Random, BenchmarkTools, SparseArrays, NonlinearEigenproblems
Random.seed!(0)
m=200;
fv=Vector{Function}(undef,m);
for i=1:m
  fv[i]=(x-> exp(i^(1/6)*x))
end
fv[1]=x->one(x); fv[2]=x->x;
Av=Vector{SparseMatrixCSC}(undef,m);
n=50;
for i=1:m
  Av[i]=sprand(n,n,0.01)
end
nep=SPMF_NEP(Av,fv);
v0=ones(n);
@btime iar(nep,maxit=100,v=v0)
