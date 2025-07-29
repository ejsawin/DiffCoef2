module SetGlobals

using CSV
using DataFrames
using Interpolations
using KernelDensity
using Statistics 


export G, M_bh, int_steps, coulomb_log, m_test, grid_size
export find_psi, find_Ftot, find_Ntot
export read_in, Ntot_grid, Ftot_grid, Ntot, Ftot, psi_itp
export m_tab, r_tab, vr_tab, vt_tab
#export id_tab, startype_tab, binflag
export vmin, vmax, rmin, rmax


#Constants
const G=1
const M_bh=8.500250825e-03
const int_steps=1000
const grid_size=0.001


# Psi via Henon

function find_psi_exact(r_, m_)
    N = length(r_)
    
    U = zeros(N + 1)
    M = zeros(N + 1)
    
    U[N + 1] = 0.0
    M[N] = sum(m_)
    
    # Potential at edge 
    U[N] = -G * (M[N] + M_bh) / r_[N]
    
    # k = N - 1 --> k = 1
    for k in (N-1):-1:1
        M[k] = M[k+1] - m_[k+1]
        U[k] = U[k+1] - G * M[k] * (1/r_[k] - 1/r_[k+1]) - G * M_bh / r_[k]
    end
    
    function psi_exact(r)
        
        k = find_shell(r, r_)

        # First shell 
        if k == 0
            return U[1] - G * M_bh * (1/r - 1/r_[1])
        
        # Last shell
        elseif k == N
            total_mass = M[N] + M_bh
            return -G * total_mass / r

        # Henon (Equation 7)
        else
            r_k = r_[k]
            r_k1 = r_[k + 1]
            U_k = U[k]
            U_k1 = U[k + 1]
            
            top = 1/r_k - 1/r
            bottom = 1/r_k - 1/r_k1
            
            return U_k + (top / bottom) * (U_k1 - U_k)
        end
    end
    
    return psi_exact
end

function find_shell(r, r_)
    N = length(r_)

    # Check edge cases 
    if r < r_[1]
        return 0
    elseif r > r_[N]
        return N
    elseif r == r_[N] 
        return N - 1
    end
    
    left = 1
    right = N - 1

    # Bisection 
    while left < right
        mid = div(left + right + 1, 2)
        if r_[mid] <= r
            left = mid
        else
            right = mid - 1
        end
    end
    
    return left
end

function find_Ntot(tab_r,tab_v,tab_m,rmin,rmax,vmin,vmax,dr,dv)
    
    Ni = Int(floor((log(rmax) - log(rmin)) / dr))
    Nj = Int(floor((log(vmax) - log(vmin)) / dv))
    
    m_sum_grid = zeros(Float64, Ni+1, Nj+1)
    
    # Bin data into cells by (i,j)
    for (r, v, m) in zip(tab_r, tab_v, tab_m)
        logr = log(r)
        logv = log(v)
        i = Int(floor((logr - log(rmin)) / dr))
        j = Int(floor((logv - log(vmin)) / dv))
        if 0 <= i <= Ni && 0 <= j <= Nj
            m_sum_grid[i+1, j+1] += m
        end
    end
    
    # Compute Ntot
    Ntot_grid = m_sum_grid ./ (dr * dv)
    
    # Create grid vectors
    r_grid = [log(rmin) + dr * i for i in 0:Ni]
    v_grid = [log(vmin) + dv * j for j in 0:Nj]
    
    function Ntot(r, v)

        #Check physicality 
        if r < rmin || r > rmax || v < vmin || v > vmax
            return 0.0
        end
        
        logr = log(r)
        logv = log(v)
        
        il = Int(floor((logr - log(rmin)) / dr))
        jl = Int(floor((logv - log(vmin)) / dv))
        
        # Clamp to valid range
        il = max(0, min(il, Ni-1))
        jl = max(0, min(jl, Nj-1))
        
        ir = il+1
        jr = jl+1
        
        ntot_ll = Ntot_grid[il+1,jl+1]
        ntot_lr = Ntot_grid[il+1,jr+1]
        ntot_rl = Ntot_grid[ir+1,jl+1]
        ntot_rr = Ntot_grid[ir+1,jr+1]
        
        # Grid coordinates
        x1 = log(rmin) + il * dr
        x2 = log(rmin) + ir * dr
        y1 = log(vmin) + jl * dv
        y2 = log(vmin) + jr * dv
        
        # Bilinear interpolation weights
        w11 = (x2 - logr) * (y2 - logv) / ((x2 - x1) * (y2 - y1))
        w12 = (x2 - logr) * (logv - y1) / ((x2 - x1) * (y2 - y1))
        w21 = (logr - x1) * (y2 - logv) / ((x2 - x1) * (y2 - y1))
        w22 = (logr - x1) * (logv - y1) / ((x2 - x1) * (y2 - y1))
        
        return w11 * ntot_ll + w12 * ntot_lr + w21 * ntot_rl + w22 * ntot_rr
    end
    
    return Ntot_grid, Ntot
end

function find_Ftot(tab_r,tab_v,tab_m,rmin,rmax,vmin,vmax,dr,dv)
    
    Ni = Int(floor((log(rmax) - log(rmin)) / dr))
    Nj = Int(floor((log(vmax) - log(vmin)) / dv))
    
    m2_sum_grid = zeros(Float64, Ni+1, Nj+1)
    
    # Bin data into cells by (i,j)
    for (r, v, m) in zip(tab_r, tab_v, tab_m)
        logr = log(r)
        logv = log(v)
        i = Int(floor((logr - log(rmin)) / dr))
        j = Int(floor((logv - log(vmin)) / dv))
        if 0 <= i <= Ni && 0 <= j <= Nj
            m2_sum_grid[i+1, j+1] += m ^ 2
        end
    end
    
    # Compute Ntot
    Ftot_grid = m2_sum_grid ./ (dr * dv)
    
    # Create grid vectors
    r_grid = [log(rmin) + dr * i for i in 0:Ni]
    v_grid = [log(vmin) + dv * j for j in 0:Nj]
    
    function Ftot(r, v)

        #Check physicality 
        if r < rmin || r > rmax || v < vmin || v > vmax
            return 0.0
        end
        
        logr = log(r)
        logv = log(v)
        
        il = Int(floor((logr - log(rmin)) / dr))
        jl = Int(floor((logv - log(vmin)) / dv))
        
        # Clamp to valid range
        il = max(0, min(il, Ni-1))
        jl = max(0, min(jl, Nj-1))
        
        ir = il+1
        jr = jl+1
        
        ftot_ll = Ftot_grid[il+1,jl+1]
        ftot_lr = Ftot_grid[il+1,jr+1]
        ftot_rl = Ftot_grid[ir+1,jl+1]
        ftot_rr = Ftot_grid[ir+1,jr+1]
        
        # Grid coordinates
        x1 = log(rmin) + il * dr
        x2 = log(rmin) + ir * dr
        y1 = log(vmin) + jl * dv
        y2 = log(vmin) + jr * dv
        
        # Bilinear interpolation weights
        w11 = (x2 - logr) * (y2 - logv) / ((x2 - x1) * (y2 - y1))
        w12 = (x2 - logr) * (logv - y1) / ((x2 - x1) * (y2 - y1))
        w21 = (logr - x1) * (y2 - logv) / ((x2 - x1) * (y2 - y1))
        w22 = (logr - x1) * (logv - y1) / ((x2 - x1) * (y2 - y1))
        
        return w11 * ftot_ll + w12 * ftot_lr + w21 * ftot_rl + w22 * ftot_rr
    end
    
    return Ftot_grid, Ftot
end

function read_in(filename)

    # read in data
    data=CSV.File(filename)
    omega_cen=DataFrame(data)
    global id_tab, m_tab, r_tab, vr_tab, vt_tab, startype_tab, binflag = eachcol(omega_cen[:,1:7])

    # set up arrays
    global v_tab = sqrt.(vr_tab .^ 2 + vt_tab .^2)
    

    # find bounds
    global vmin=minimum(v_tab)
    global vmax=maximum(v_tab)

    global rmin=minimum(r_tab)
    global rmax=maximum(r_tab)

    global coulomb_log = M_bh / mean(m_tab)
    global m_test = mean(m_tab) # Test mass

    #find df, potential 
    global psi_itp=find_psi_exact(r_tab,m_tab)
    global Ntot_grid,Ntot=find_Ntot(r_tab,v_tab,m_tab,rmin,rmax,vmin,vmax,grid_size,grid_size)
    global Ftot_grid,Ftot=find_Ftot(r_tab,v_tab,m_tab,rmin,rmax,vmin,vmax,grid_size,grid_size)
end

end