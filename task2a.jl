
using Plots
include("directdiag.jl")

function print_state(s)
    res  = "|"
    for c in s
        if c == 0
            res *= "↓"
        elseif c == 1
            res *= "↑"
        end
    end
    res *= "⟩"
    println(res)
end

function print_state(s)
    res  = "|"
    for c in s
        if c == 0
            res *= "$c"
        elseif c == 1
            res *= "$c"
        end
        res *= ";"
    end
    res *= "⟩"
    println(res)
end

function task()
    E0 = []
    eigenvecs = []
    mzs = []


    for N in 2:6
    
        H = createHamiltonian(N)
        res = eigen(H)
        push!(E0, res.values[1])
        push!(eigenvecs, res.vectors[:,1])
        # push!(mzs, sum([res.vectors[i,1] * [Sᶻ(BittArray(digits(j, 2, pad=N)), j) for j in 1:N] for i in 1:N]))
        println(size(res.vectors))
    end
    return E0, eigenvecs, mzs
end

res = task()

pltEo = scatter(res[1])

pltmzs = scatter(res[3])

for i in 1:5
    print_state(res[2][i])
end

# display(pltEo)
# display(pltmzs)
println(res)


