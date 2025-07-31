# Orbit bounds using manual bisection 
function orbit_interp(E_, L_)
    f(x) = psi_itp(x) + (L_^2)/(2*x^2) - E_
    
    # Find inner turning point
    r_inner = bisection(f,rmin,(rmin+rmax)/2,tol,bisection_steps)
    
    # Find outer turning point  
    r_outer = bisection(f,(rmin+rmax)/2,rmax,tol,bisection_steps)
    
    return [r_inner, r_outer]
end

