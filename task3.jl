include("directdiag.jl")

using LinearAlgebra
using LaTeXStrings
using Plots

Ns = 2:13

function task(N)
    bloxs = allmzBlocks(N)
    Emin = Inf
    for blok in bloxs
        E0 = eigen(blok).values[1]
        Emin = min(E0,Emin)
    end
    return Emin
end

x_old=[@elapsed(eigen(createHamiltonian(n))) for n in Ns]
x=[@elapsed(task(n)) for n in Ns]

plt = scatter(Ns,x_old,
    xaxis=L"Number of spins $N$",
    yticks= [0.0001,0.001,0.01,0.1,1,10,100, 1000],
    xticks= 1:14,
    yaxis="Time for diagonalisation [s]",
    yscale=:log10,
    label="Full Hamiltonian",
    legend=:bottomright
)



scatter!(Ns,x,
    xaxis=L"Number of spins $N$",
    yticks= [0.0001,0.001,0.01,0.1,1,10,100, 1000],
    xticks= 1:14,
    yaxis="Time for diagonalisation [s]",
    yscale=:log10,
    label=L"With $m_z$ blocks"
)


savefig(plt, "saves/task3.timedependenc.pdf")
display(plt)


# Last points: 8sec vs 40 sec.