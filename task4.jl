include("directdiag.jl")
using Plots
using ProgressMeter

function task(N)
    blocks = allmzBlocks(N)
    #H = createHamiltonian(N)

    E0 = Inf
    @showprogress for block in blocks
        M = size(block)[1]
        res = eigen(lanzos(block, M))
        E0 = min(E0, res.values[1])
    end
    return E0
end


task(5)
