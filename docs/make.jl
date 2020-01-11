using Documenter, MolecularConformation

makedocs(sitename= "MolecularConformation.jl Documentation",
				   pages = ["Home"=> "index.md",
							"Quick Start" => "quick.md",
							"API" => "api.md",
							"Contributing"=> "contributing.md",
							"Development Status"=> "devstatus.md"],
				   format = Documenter.HTML(prettyurls = false)
				   
		)


