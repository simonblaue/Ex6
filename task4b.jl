include("directdiag.jl")


function task(Λ)
    
    mz_blocks = allmzBlocks(10)
    mz0_block = mz_blocks[floor(Int,length(mz_blocks)/2)]

    H = lanzos(mz0_block, Λ)
    
    return eigen(H).values[1:4]
end

task.(4:30)