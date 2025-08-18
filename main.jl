# Necessary modules
include("make_globals.jl")
include("loc_vel_coef.jl")
include("vel_calc.jl")
include("midpoint.jl")

using .MakeGlobals
using Plots
using KernelDensity


# Read in data
read_file("plummer_IC.csv")


# Test
r_=1.5
vr=0.7
vt=0.3

E_test=0.5(vr^2+vt^2)+psi_calc(r_)
L_test=vt/r_

println(find_locv(r_,E_test,L_test))



# Graphing

#r_test=exp.(range(log(rmin),log(rmax),4000))
#v_test=exp.(range(log(vmin),log(vmax),4000))

#kde1=kde((r_tab,v_tab))

#graph=heatmap(r_test,v_test,F_1D,xlim=(0,4),title="Interpolated DF")
#graph2=heatmap(r_test,v_test,N_1D,xlim=(0,4),title="Interpolated DF")
#graph3=plot(graph,graph2,layout=(1,2))
#display(graph3)
#readline()