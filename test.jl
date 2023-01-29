
include("directdiag.jl")

@time H = createHamiltonian(2)

res = eigen(H)
