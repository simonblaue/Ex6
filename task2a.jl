
using Plots
include("directdiag.jl")


function task()
    E0 = []
    eigenvecs = []
    mzs = []


    for N in 2:6
    
        H = createHamiltonian(N)
        res = eigen(H)
        push!(E0, res.values[1])
        push!(eigenvecs, res.vectors[:,1])
        #push!(mzs, sum([res.vectors[i,1] * [Sᶻ(BittArray(digits(j, 2, pad=N)), j) for j in 1:N] for i in 1:N]))
    end
    return E0, eigenvecs, mzs
end

res = task()

pltEo = scatter(res[1])

pltmzs = scatter(res[3])


display(pltEo)
# display(pltmzs)



