using LinearAlgebra
using BitOperations

function createHamiltonian(N)
    H = zeros(2^N,2^N)
    J = 1
    for a in 0:2^N-1
        for i in 0:N-1
            j = mod(i+1,N)
            if bget(a,i) == bget(a,j)
                H[a+1,a+1] += 1/4
            else
                H[a+1,a+1] += -1/4
                b = bflip(a, i)
                b = bflip(b, j)
                H[a+1, b+1] += 1/2
            end 
        end
    end
    H .*= J 
    return H
end

function m_z(vector::Vector{Float64}, N::Int)
    base = [BitArray(digits(i-1, base=2, pad=N)) for i in 1:2^N]
    S(state) = [s ? 1/2 : -1/2 for s in state]
    mz = sum(sum(vector .* S.(base)))
    return mz
end

### Task 3

function sa(N, nup)
    sₐ = []
    for s in 0:2^N-1
        if sum(BitArray(digits(s, base=2, pad=N))) == nup
            push!(sₐ,s)
        end
    end
    return sₐ
end

function createBlockHamiltonian(states)
    H_size = length(states)
    N = length(BitArray(digits(maximum(states), base=2)))
    H = zeros(H_size,H_size)
    for (k,a) in enumerate(states)
        for i in 0:N-1
            j = mod(i+1,N)
            if bget(a,i) == bget(a,j)
                H[k,k] += 1/4
            else
                H[k,k] += -1/4
                aflipped = bflip(a, i)
                aflipped = bflip(aflipped, j)
                l = findall(aflipped .== states)[1]
                H[k, l] += 1/2
            end 
        end
    end
    return H
end

function allmzBlocks(N)
    return [createBlockHamiltonian(sa(N,nup)) for nup in 0:N]
end

#### Task 4

function lanzos(H, Λ)
    
    @assert size(H)[1] == size(H)[2] "Not a square Matrix"

    N = size(H)[1]
    
    if Λ == -1
        Λ = N
    end

    if N < 2
        return H
    end

    f = zeros(N, Λ)
    a = ones(Λ)
    b = ones(Λ-1)

    f[:,1] = rand(N)

    a[1] = (f[:,1]' * H * f[:,1])/ (f[:,1]⋅f[:,1])

    f[:,2] = H * f[:,1] - a[1] * f[:,1]

    λ = 2
    while λ < Λ

        a[λ] = (f[:,λ]' * H *f[:,λ]) / (f[:,λ]⋅f[:,λ])
        b[λ-1] = (f[:,λ]⋅f[:,λ]) / (f[:,λ-1]⋅f[:,λ-1])

        f[:,λ+1] = H * f[:,λ] - a[λ] * f[:,λ] - b[λ-1] * f[:,λ-1]

        λ += 1
    end

    smallH = SymTridiagonal(a, sqrt.(b))

    return smallH
end