include("directdiag.jl")

using LinearAlgebra
using LaTeXStrings
using Plots

Ns = 2:11

function task(N)
    bloxs = [createBlockHamiltonian(sa(N,nup), N) for nup in 0:N]
    Emin = Inf
    for blok in bloxs
        if size(blok) != 0
            E0 = eigen(blok).values[1]
            Emin = min(E0,Emin)
        end
    end
    return Emin
end


x=[@elapsed(task(n)) for n in Ns]

plt = scatter(Ns,x,
    xaxis=L"Number of spins $N$",
    yticks= [0.0001,0.001,0.01,0.1,1,10,100, 1000],
    xticks= 1:14,
    yaxis="Time for diagonalisation [s]",
    yscale=:log10,
    legend=false
)


savefig(plt, "saves/task3.timedependenc.pdf")
display(plt)