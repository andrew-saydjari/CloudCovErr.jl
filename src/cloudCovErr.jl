module cloudCovErr

include("cov.jl")

include("perstar.jl")

include("preprocess.jl")

include("decam.jl")
using .decam

#include("valid.jl")

# include("plotting.jl")
# using .plotting

end
