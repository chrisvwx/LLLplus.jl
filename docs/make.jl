using LLLplus, Documenter, MacroTools

makedocs(
         sitename="LLLplus.jl",
         pages = [
                  "Home" => "index.md",
                  "Functions" => "functions.md"],
         format = Documenter.HTML(prettyurls = haskey(ENV, "CI")))

