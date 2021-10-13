using disCovErr
using Documenter

makedocs(
    clean=true,
    highlightsig = true,
    sitename= "disCovErr",
)

deploydocs(
    repo = "github.com/andrew-saydjari/disCovErr.jl.git",
    branch = "gh-pages",
)
