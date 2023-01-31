
include("directdiag.jl")

# N = 6

# @time H = createHamiltonian(N)

# res = eigen(H)

# vector = res.vectors[:,1]
# # display(vector)

# function correlation(eigenvec)
#     base = [BitArray(digits(i-1, base=2, pad=N)) for i in 1:2^N]
#     S(state) = [s ? 1/2 : -1/2 for s in state]
#     Sz = vector .* S.(base)
#     corr = 0
#     for bv in Sz
#         corr += bv[1]*bv[3]
#     end
#     return corr
# end

# sum(correlation.(res.vectors))

N = 1

mod1(3,N)