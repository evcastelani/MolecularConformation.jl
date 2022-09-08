using Plots

#pgfplotsx()

function perfomance_profile(c, method = ["Alg $(i)" for i=1:size(c,2)], color = [:black,:gray],style = [:dot,:solid])
    cmin = minimum(c,dims=2)
    R = c ./ cmin
    t = sort(unique(R))
    if t[end] == NaN
        pop!(t)
    end
    plot(title = " ", xlabel = "Perfomance Ratio", ylabel = "Problems solved")
    plot!(xticks=round.(0.0:0.5:3,digits=2))
    plot!(xaxis=:log)
    for i = 1:size(c, 2)
        plot!(t, [sum(R[:,i] .<= ti)/size(c,1) for ti in t], label=method[i], t=:steppre, lw=2,linecolor = color[i],line = style[i])
    end
    ylims!(0, 1)
    #plot!(xformatter = xi -> "\$10^{$xi}\$")
    savefig("tex/figures/tp_perf.pdf")
end


function filtering(c)
    r = zeros(size(c))
    k = 1
    for i=1:length(c[:,1])
        if c[i,1] == Inf && c[i,2] != Inf
            println("problem $(i) were solved by 2 but does not 1")
        elseif c[i,1] != Inf && c[i,2] == Inf
            println("problem $(i) were solved by 1 but does not 2")
        elseif c[i,1] == Inf && c[i,2] ==Inf
            println("problem $(i) were not solved! This line was deleted")
        else
            r[k,:] = c[i,:]
            k = k+1
        end
    end
    return Float64.(r[1:k-1,:])
end