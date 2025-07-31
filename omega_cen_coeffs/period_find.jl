# Period calculation via anomaly 

function period_integrand(tt, _E, _L, rp, ra)
    r_hat = (ra - rp)/2
    t_hat = (ra + rp)/2
    r = r_hat * sin(tt) + t_hat
    vr = v_r(r, _E, _L)
    return r_hat * cos(tt) / vr
end

function period(_E, _L)
    rp, ra = orbit_interp(_E, _L)
    return 2 * midpoint(tt -> period_integrand(tt, _E, _L, rp, ra), -pi/2, pi/2, int_steps)
end