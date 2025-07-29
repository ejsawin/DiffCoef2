include("set_globals.jl")
include("midpoint.jl")
include("velocities.jl")
include("local_coef.jl")
include("orbits.jl")
include("period_find.jl")
include("local_iom.jl")
include("average_coef.jl")
include("visualizer.jl")

using Plots
using LaTeXStrings
using QuadGK
using KernelDensity
using Cubature
using Roots
using .SetGlobals

SetGlobals.read_in("omegaCenEddie.csv")

#r_test=range(rmin,rmax,1000) # Test arrays
#v_test=range(vmin,vmax,1000)

#res = [log10.(Ntot(y,z)) for z in v_test, y in r_test] # Matrix 
#res2 = [log10.(Ftot(w,o)) for o in v_test, w in r_test]
#res3 = [log10.(Ftot(w,o)/Ntot(w,o)) for o in v_test, w in r_test]

#graph1 = heatmap(r_test,v_test,res,color=:tofino,size=(1500,750),title="Ntot",xlim=(0,15),ylim=(0,1.5)) # Plotting 

average_coef_plot(-0.2,0.1,1.5,100)