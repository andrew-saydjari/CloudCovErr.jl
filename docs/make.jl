using disCovErr
using Documenter

makedocs(
    root= "/",
    source= "src",
    build="build",
    clean=true,
    highlightsig = true,
    sitename= "disCovErr",
)

deploydocs(
    repo = "github.com/andrew-saydjari/disCovErr.jl.git",
)
