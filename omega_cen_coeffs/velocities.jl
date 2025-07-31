function v_t(r,E,L)
    var=L/r
    return var
end

function v_r(r,E,L)
    var1=sqrt(abs(2*(E-psi_itp(r))-L^2/r^2))
    return var1
end