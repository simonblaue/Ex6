include("directdiag.jl")
using Plots

function task4(N)
    blocks = allmzBlocks(N)
    E0 = Inf
    for block in blocks
        M = size(block)[1]
        res = eigen(Matrix(lanzos(block, M+1)))
        E0 = min(E0, res.values[1])
    end
    return E0
end


function task3(N)
    bloxs = allmzBlocks(N)
    Emin = Inf
    for blok in bloxs
        E0 = eigen(blok).values[1]
        Emin = min(E0,Emin)
    end
    return Emin
end

######################

Ns = 2:12

@time x2=[@elapsed(eigen(createHamiltonian(n))) for n in Ns]
@time x3=[@elapsed(task3(n)) for n in Ns]
@time x4=[@elapsed(task4(n)) for n in Ns]

plt = scatter(Ns,x2,
    xaxis=L"Number of spins $N$",
    yticks= [0.0001,0.001,0.01,0.1,1,10,100, 1000],
    xticks= 1:14,
    yaxis="Time for diagonalisation [s]",
    yscale=:log10,
    label="Full Hamiltonian",
    legend=:bottomright
)



scatter!(Ns,x3,
    xaxis=L"Number of spins $N$",
    yticks= [0.0001,0.001,0.01,0.1,1,10,100, 1000],
    xticks= 1:14,
    yaxis="Time for diagonalisation [s]",
    yscale=:log10,
    label=L"With $m_z$ blocks"
)

scatter!(Ns,x4,
    xaxis=L"Number of spins $N$",
    yticks= [0.0001,0.001,0.01,0.1,1,10,100, 1000],
    xticks= 1:14,
    yaxis="Time for diagonalisation [s]",
    yscale=:log10,
    label=L"With $m_z$ blocks and Lanzos"
)


savefig(plt, "saves/task4.timedependenc.pdf")
display(plt)
