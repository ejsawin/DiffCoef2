function orbit_interp(E_,L_)    
    r1=rmin
    r2=rmax

    f(x)=psi_itp(x)+(L_^2)/(2*x^2) - E_ 
    rs=find_zeros(f,r1,r2) #Find turning points 
    return rs
end