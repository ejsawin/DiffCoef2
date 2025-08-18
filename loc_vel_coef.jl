# Calculate local velocity diffusion coefficients at (r, E, L)

function find_locv(rr,EE,LL)

    # determine velocity
    v_step=sqrt((v_r(rr,EE,LL))^2 + (v_t(rr,EE,LL))^2)

    # constants
    c1=-(G^2*coulomb_log)/(v_step^2 * rr^2)
    c2=(2*G^2*coulomb_log)/(3*rr^2)

    # integrands
    vp_int(vv)=F_1D(rr,vv)
    vp_int1(vv)=N_1D(rr,vv)

    vpsq_int(vv)=(vv^2)/(v_step^3)*F_1D(rr,vv)
    vpsq_int1(vv)=(1 / vv) * F_1D(rr,vv)

    # calculate integrals
    vp_calc=midpoint(x->vp_int(x),vmin,v_step,int_steps)
    vp_calc1=midpoint(x->vp_int1(x),vmin,v_step,int_steps)

    vpsq_calc=midpoint(x->vpsq_int(x),vmin,v_step,int_steps)
    vpsq_calc1=midpoint(x->vpsq_int1(x),v_step,vmax,int_steps)

    # calculate diffusion coefficients 
    D_par=c1*(vp_calc+m_test*vp_calc1)
    D_par2=c2*(vpsq_calc+vpsq_calc1)
    D_tan2=c2*((3/v_step^2)*vp_calc - vpsq_calc + 2 * vpsq_calc1)

    return D_par,D_par2,D_tan2
end

