using disCovErr
using Documenter

makedocs(
    root= "/Users/saydjari/Dropbox/GradSchool_AKS/Doug/myPublicPkgs/disCovErr.jl/docs",
    source= "/Users/saydjari/Dropbox/GradSchool_AKS/Doug/myPublicPkgs/disCovErr.jl/src",
    build="build",
    clean=true,
    highlightsig = true,
    sitename= "disCovErr",
)

deploydocs(
    repo = "github.com/andrew-saydjari/disCovErr.jl.git",
)
