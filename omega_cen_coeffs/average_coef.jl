function average_coef_int(tt,EE,LL,ii,rp,ra)
    r_hat=(ra-rp)/2
    t_hat=(ra+rp)/2

    rrr=r_hat*sin(tt)+t_hat
    vr_tilde=v_r(rrr,EE,LL)/(r_hat*cos(tt))

    return (local_iom(rrr,EE,LL)[ii])/vr_tilde
end

function average_coef(EE_,LL_,indx)
    rp,ra=orbit_interp(EE_,LL_)
    T=period(EE_,LL_)
    return (2/T)*midpoint(x -> average_coef_int(x,EE_,LL_,indx,rp,ra),-pi/2,pi/2,int_steps)
end