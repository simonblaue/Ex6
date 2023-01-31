
include("directdiag.jl")
using Plots
using LaTeXStrings
using ProgressMeter

function correlations(eigenvec,eigenval, N)
    base = [BitArray(digits(i-1, base=2, pad=N)) for i in 1:2^N]
    S(state) = [s ? 1/2 : -1/2 for s in state]
    Sz = eigenvec .* S.(base) .* exp(-eigenval)
    corr = 0
    for bv in Sz
        corr += bv[mod1(1,N)]*bv[mod1(3,N)]
    end
    return corr
end

function averages(N)

    H = createHamiltonian(N)
    res = eigen(H)

    λ = res.values
    v = [res.vectors[:,i] for i in 1:2^N]

    Z = sum(exp.(-λ))

    cor = sum(correlations.(v,λ,N)) / Z
    m = sum(m_z.(v, N) .* exp.(-λ) ) / Z
    m² = sum(m_z.(v,N).^2 .* exp.(-λ)) / Z
    E = sum(λ .* exp.(-λ)) / Z
    E² = sum(λ.^2 .* exp.(-λ)) / Z

    Χ = m² - m^2
    C = E² - E^2

    return cor, Χ, C
end

function task(Ns)

    corrs, Χs, Cs = [], [], []

    @showprogress for N in Ns
        res = averages(N)
        push!(corrs, res[1])
        push!(Χs, res[2])
        push!(Cs, res[3])
    end
    return corrs, Χs, Cs
end

############

Ns = 4:12

corrs, Χs, Cs = task(Ns)

plt1 = plot(Ns, corrs,
    xlabel=L"Number of spins $N$",
    ylabel="Correlation",
    label=""
)

plt2 = plot(Ns, Χs,
    xlabel=L"Number of spins $N$",
    ylabel="Magnetic susceptibility",
    label=""
)
hline!([1/4], color="black", label="high-temperature limit")

plt3 = plot(Ns, Cs,
    xlabel=L"Number of spins $N$",
    ylabel="specific heat",
    label=""
)
hline!([3/13], color="black", label="high-temperature limit")

display(plt1)
display(plt2)
display(plt3)

savefig(plt1, "saves/task2c.corr.pdf")
savefig(plt2, "saves/task2c.suzept.pdf")
savefig(plt3, "saves/task2c.specheat.pdf")
