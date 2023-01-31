using LinearAlgebra
using BitOperations

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

