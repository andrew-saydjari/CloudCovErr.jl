using cloudCovErr
using Documenter

makedocs(
    modules = [cloudCovErr, decam],
    clean=true,
    highlightsig = true,
    sitename= "cloudCovErr.jl",
    format = Documenter.HTML(
        assets = ["assets/favicon.ico"],
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    pages    = [
        "Introduction" => "index.md",
        "API Reference" => "api.md",
        "Contributing" => "contrib.md"
    ]
)

deploydocs(
    repo = "github.com/andrew-saydjari/cloudCovErr.jl.git",
    branch = "gh-pages",
    devbranch = "main"
)
