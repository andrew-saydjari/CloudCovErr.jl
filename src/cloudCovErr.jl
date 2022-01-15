module CloudCovErr

include("cov.jl")

include("perstar.jl")

include("preprocess.jl")

include("decam.jl")
using .decam

end
