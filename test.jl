include("directdiag.jl")


Ns = 2:6

function task(N)
    bloxs = allmzBlocks(N)
    Emin = Inf
    for blok in bloxs
        E0 = eigen(blok).values[1]
        Emin = min(E0,Emin)
    end
    return Emin
end


# task(6)

# bloxs = size.(allmzBlocks(6))

N = 6

states_all = [sa(N,nup) for nup in 0:N]

Nsmall = [createBlockHamiltonian(states) for states in states_all]
