using Documenter, MolecularConformation

makedocs(sitename= "MolecularConformation.jl",
				   pages = ["Home"=> "index.md",
							"Quick Start" => "quick.md",
							"API" => "api.md",
							"Contributing"=> "contributing.md",
							"Development Status"=> "devstatus.md"],
				   format = Documenter.HTML(prettyurls = false)
				   
		)

deploydocs(
		   	repo = "github.com/evcastelani/MolecularConformation.jl.git"
				)
