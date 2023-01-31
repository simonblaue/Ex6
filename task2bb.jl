include("directdiag.jl")

using LinearAlgebra
using LaTeXStrings
using Plots

Ns = 2:13

x=[@elapsed(eigen(createHamiltonian(n))) for n in Ns]

plt = scatter(Ns,x,
    xaxis=L"Number of spins $N$",
    yticks= [0.0001,0.001,0.01,0.1,1,10,100, 1000],
    xticks= 1:14,
    yaxis="Time for diagonalisation [s]",
    yscale=:log10,
    legend=false
)


savefig(plt, "saves/task2b.timedependenc.pdf")
display(plt)