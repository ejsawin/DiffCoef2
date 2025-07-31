function local_coef(rr,EE,LL)

    # Find step velocity, constants
    v_step=sqrt((v_r(rr,EE,LL)) ^ 2 + (v_t(rr,EE,LL)) ^ 2)

    c1=-((G^2)*coulomb_log)/((v_step^2)*(rr^3))
    c2=(2*(G^2)*coulomb_log)/(3*(rr^3))

    # Define, evaluate integrals, everything in log spacing. 
    vp_int(vv)=Ntot_log_1d(rr,vv)
    vp_int1(vv)=Ftot_log_1d(rr,vv)

    vpsq_int(vv)=((vv^2)/(v_step^3))*Ftot_log_1d(rr,vv)
    vpsq_int1(vv)=(1/(vv))*Ftot_log_1d(rr,vv)

    vp_calc=midpoint(x->vp_int(x),vmin,v_step,int_steps)
    vp1_calc=midpoint(x->vp_int1(x),vmin,v_step,int_steps)

    vpsq_calc=midpoint(x->vpsq_int(x),vmin,v_step,int_steps)
    vpsq1_calc=midpoint(x->vpsq_int1(x),v_step,vmax,int_steps)

    # Evaluate coefficients 
    D_par=c1*(vp_calc+m_test*vp1_calc)
    D_parsq=c2*(vpsq_calc+vpsq1_calc)
    D_tansq=c2*((3/v_step)*vp1_calc - vpsq_calc + 2*vpsq1_calc)

    return D_par, D_parsq, D_tansq
end
