function local_coef(rr,EE,LL)

    # Calculate step velocity
    v_step=sqrt((v_r(rr,EE,LL)) ^ 2 + (v_t(rr,EE,LL)) ^ 2)

    c1 = -16*(pi^2)*(G^2)*coulomb_log/(v_step^2)
    c2 = 32*(pi^2)*(G^2)*coulomb_log/3

    vp_int(vv) = (1/(16*(pi^2)*(vv^2)*(rr^2)))*(vv^2)*(Ftot(rr,vv) + m_test * Ntot(rr,vv))

    vp2_int(vv) = (1/(16*(pi^2)*(vv^2)*(rr^2)))*(vv^4)*(1/v_step^3)*Ftot(rr,vv)
    
    vp2_int1(vv) = (1/(16*(pi^2)*(vv^2)*(rr^2)))*vv*Ftot(rr,vv)

    vt2_int(vv) = (1/(16*(pi^2)*(vv^2)*(rr^2)))*((3*vv^2)*(1/v_step) - (vv^4)*(1/v_step^3))*Ftot(rr,vv)
    
    vt2_int1(vv) = (1/(16*(pi^2)*(vv^2)*(rr^2)))*2*vv*Ftot(rr,vv)

    vp=c1*midpoint(x -> vp_int(x),vmin,v_step,int_steps)
    vp2=c2*(midpoint(x -> vp2_int(x),vmin,v_step,int_steps)+midpoint(x -> vp2_int1(x),v_step,vmax,int_steps))
    vt2=c2*(midpoint(x -> vt2_int(x),vmin,v_step,int_steps)+midpoint(x -> vt2_int1(x),v_step,vmax,int_steps))

    return vp,vp2,vt2
end
