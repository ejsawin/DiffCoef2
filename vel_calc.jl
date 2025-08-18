# Calculate velocities for (r,E,L)

# Radial velocity- Eq 1.3
function v_r(rr,EE,LL)
    return sqrt(abs(2*(EE-psi_calc(rr))-(LL^2)/(rr^2)))
end

# Tangential velocity - Eq 1.2
function v_t(rr,EE,LL)
    return LL/rr
end