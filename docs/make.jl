using cloudCovErr
using Documenter

makedocs(
    clean=true,
    highlightsig = true,
    sitename= "cloudCovErr.jl",
    pages    = [
        "Introduction" => "index.md",
        "API Reference" => "api.md"
    ]
)

deploydocs(
    repo = "github.com/andrew-saydjari/cloudCovErr.jl.git",
    branch = "gh-pages",
    devbranch = "main"
)
