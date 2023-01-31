
include("directdiag.jl")

N = 2

@time H = createHamiltonian(N)

res = eigen(H)

base = [BitArray(digits(i-1, base=2, pad=N)) for i in 1:2^N]

S(state) = [s ? 1/2 : -1/2 for s in state]

(sum(res.vectors[1] * S.(base)))

H
