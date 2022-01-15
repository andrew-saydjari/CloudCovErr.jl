using CloudCovErr
using Documenter

makedocs(
    modules = [CloudCovErr, CloudCovErr.decam],
    clean=true,
    highlightsig = true,
    sitename= "CloudCovErr.jl",
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
    repo = "github.com/andrew-saydjari/CloudCovErr.jl.git",
    branch = "gh-pages",
    devbranch = "main"
)
