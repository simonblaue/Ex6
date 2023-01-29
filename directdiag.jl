using LinearAlgebra
using BitOperations

function ⋆(s1,s2)
    if s1 == 0 || s2 == 0
        return 0
    else
        return dot(s1,s2)
    end
end


function S⁺(s::Int, i)
    return s
end

function S⁻(s::Int, i)
    return s
end


function S⁺(s::BitArray, i)
    if s[i] == 1
        return 0
    else
        s_new = copy(s)
        s_new[i] = 1
        return s_new
    end
end


function S⁻(s::BitArray, i)
    if s[i] == 0
        return 0
    else
        s_new = copy(s)
        s_new[i] = 0
        return s_new
    end
end

function Sᶻ(s, i)
    return s[i]==1 ? 1/2 : -1/2
end

function mz(s)
    
end

function createHamiltonian(N, s::String)
    H = zeros(2^N,2^N)
    J = 1
    for n in 1:2^N
        for m in 1:2^N
            s1 = BitArray(digits(n-1, base=2, pad=N))
            s2 = BitArray(digits(m-1, base=2, pad=N))
            for i in 1:N
                j = mod1(i+1,N)
                H[n,m] += Sᶻ(s2,i) * Sᶻ(s2,j) * (s1 ⋆ s2)
                H[n,m] += 1/2 * (s1 ⋆ S⁺(S⁻(s2, j), i))
                H[n,m] += 1/2 * (s1 ⋆ S⁻(S⁺(s2, j), i))
            end
            H[n,m] *= J 
        end
    end
    return H
end


function createHamiltonian(N)
    H = zeros(2^N,2^N)
    J = 1
    for a in 0:2^N-1
        for i in 0:N-1
            j = mod(i+1,N)
            if bget(a,i) == bget(a,j)
                H[a+1,a+1] = 1/4
            else
                H[a+1,a+1] = -1/4
                b = bflip(a, i)
                b = bflip(b, j)
                H[a+1, b+1] = 1/2
            end 
        end
    end
    H .*= J 
    return H
end



