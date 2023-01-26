
include("directdiag.jl")

H = createHamiltonian(2)

res = eigen(H)

H * res.vectors[:,1]