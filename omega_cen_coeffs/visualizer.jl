function average_coef_plot(E_fixed,lower,upper,steps)
    ll_range=exp.(range(log(lower),log(upper),length=steps))

    de=zeros(Float64,length(ll_range))
    de2=zeros(Float64,length(ll_range))
    dl=zeros(Float64,length(ll_range))
    dl2=zeros(Float64,length(ll_range))
    del=zeros(Float64,length(ll_range))

    Threads.@threads for i in eachindex(ll_range)
        de[i] = abs(average_coef(E_fixed,ll_range[i],1))
        de2[i] = abs(average_coef(E_fixed,ll_range[i],2))
        dl[i] = abs(average_coef(E_fixed,ll_range[i],3))
        dl2[i] = abs(average_coef(E_fixed,ll_range[i],4))
        del[i] = abs(average_coef(E_fixed,ll_range[i],5))
    end

    graph1=plot(ll_range,[de de2 dl dl2 del],label=["DE" "DE2" "DL" "DL2" "DEL"],title="Omega Cen Diff Coefs",xscale=:log10,yscale=:log10,xlabel="log_10(L)",ylabel="log_10")
    display(graph1)
    savefig("coef_plot.png")
    readline()
end