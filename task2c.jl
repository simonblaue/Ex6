
include("directdiag.jl")


function correlation_time_v(eigenvec)
    base = [BitArray(digits(i-1, base=2, pad=N)) for i in 1:2^N]
    S(state) = [s ? 1/2 : -1/2 for s in state]
    Sz = vector .* S.(base)
    corr = 0
    for bv in Sz
        corr += bv[1]*bv[3]
    end
    return corr
end

function averages(N)

    H = createHamiltonian(N)
    res = eigen(H)

    λ = res.values
    v = res.vectors

    Z = sum(exp.(-λ))





    return Z
end

N = 2

averages(N)