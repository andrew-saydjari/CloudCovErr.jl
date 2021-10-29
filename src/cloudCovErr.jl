module cloudCovErr

include("cov.jl")

include("perstar.jl")

include("preprocess.jl")

include("decam.jl")
using .decam

# include("plotting.jl")
# using .plotting

end
