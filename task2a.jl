
using Plots
using Printf
using LaTeXStrings
include("directdiag.jl")

function task(Ns)
    E0 = []
    eigenvecs = []
    mzs = []


    for N in Ns
    
        H = createHamiltonian(N)
        res = eigen(H)
        push!(E0, res.values[1])
        push!(eigenvecs, res.vectors[:,1])
        push!(mzs, m_z(res.vectors[:,1], N))
    end
    return E0, eigenvecs, mzs
end

function eigvec_LaTeX(vec)
    res = "\$("

    for num in vec
        if num == 0
            res *= "0, "
        else
            res *= "\\SI{"
            res *= @sprintf "%.2e" num
            res *= "}{}, "
        end
    end

    res *= ")^\\top\$"
    return res
end


###########

Ns = 2:6

res = task(Ns)

pltEo = scatter(Ns,res[1],
    xlabel = L"Number of Spins $N$",
    ylabel = L"Ground-state Energy $E_0$",
    label=""
)

pltmzs = scatter(Ns,res[3],
    xlabel = L"Number of Spins $N$",
    ylabel = L"Magnetisation $m_z$",
    label=""
)


savefig(pltEo, "saves/task2a.E0.pdf")
savefig(pltmzs, "saves/task2a.mzs.pdf")

display(pltEo)
display(pltmzs)

for v in res[2]
    println(eigvec_LaTeX(v) * "\\\\")
    println()
end




