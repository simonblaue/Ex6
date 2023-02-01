include("directdiag.jl")


@time blocks = allmzBlocks(4)

size.(blocks)