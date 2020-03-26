using LLLplus, Documenter, MacroTools

setup = quote
    using LLLplus
end
DocMeta.setdocmeta!(LLLplus, :DocTestSetup, setup; recursive = true)

makedocs(
         sitename="LLLplus.jl",
         pages = [
                  "Home" => "index.md",
                  "Functions" => "functions.md"],
         format = Documenter.HTML(prettyurls = haskey(ENV, "CI")),
         strict = true,
)

