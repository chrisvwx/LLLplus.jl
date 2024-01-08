using LLLplus, Documenter, MacroTools

# make.jl still doesn't seem to be running doctests.  Or at least it isn't
# running the doctests I expect it to run
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
)

deploydocs(
    repo   = "github.com/chrisvwx/LLLplus.jl.git"
)
