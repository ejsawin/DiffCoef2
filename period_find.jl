# Still need to convert to anomaly! 

function period(_E,_L)
    rp_calc,ra_calc=orbit_interp(_E,_L)
    T_calc=2*midpoint(x -> 1/v_r(x,_E,_L),rp_calc,ra_calc,int_steps)
    return T_calc
end