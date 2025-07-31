module SetGlobals

using CSV
using DataFrames
using Interpolations
using KernelDensity
using Statistics 


export G, M_bh, int_steps, coulomb_log, m_test, grid_size, bisection_steps, tol
export m_tab, r_tab, vr_tab, vt_tab
export vmin, vmax, rmin, rmax
export find_psi_arrays, find_psi_exact, psi_exact, psi_itp, U_tab, M_tab
export read_in, Ntot_grid, Ftot_grid, Ntot_log_1d, Ftot_log_1d

# Constants
const G=1
const M_bh=0.0 #8.500250825e-03

# Numerical 
const int_steps=1000 # Integration steps 
const grid_size=0.01 # Grid size, log

const bisection_steps=10 #Bisection parameters
const tol=1e-8

# Potential Calculation

# --- Shell finding via Bisection ---
function find_shell(r, r_)
    N = length(r_)

    # Check bounds
    if r < r_[1]
        return 0
    elseif r > r_[N]
        return N
    elseif r == r_[N]
        return N - 1
    end

    # Bisection 
    left, right = 1, N - 1
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

#Generate U,M arrays 
function find_psi_arrays(r_, m_)
    N = length(r_)
    U = zeros(N + 1)
    M = zeros(N + 1)

    U[N + 1] = 0.0
    M[N] = sum(m_)
    U[N] = -G * (M[N] + M_bh) / r_[N]

    for k in (N - 1):-1:1
        M[k] = M[k + 1] - m_[k + 1]
        U[k] = U[k + 1] - G * M[k] * (1/r_[k] - 1/r_[k + 1]) - G * M_bh / r_[k]
    end

    return U, M
end

# --- Compute psi via  Henon ---
function find_psi_exact(r, r_, U, M)
    N = length(r_)
    k = find_shell(r, r_)

    # First shell
    if k == 0
        return U[1] - G * M_bh * (1/r - 1/r_[1])

    # Last shell 
    elseif k == N
        return -G * (M[N] + M_bh) / r

    # Between 
    else
        rk, rk1 = r_[k], r_[k + 1]
        Uk, Uk1 = U[k], U[k + 1]
        return Uk + ((1/rk - 1/r) / (1/rk - 1/rk1)) * (Uk1 - Uk)
    end
end

# --- Wrapper to make callable psi(r)---
function psi_exact(U, M, r_)
    return r -> find_psi_exact(r, r_, U, M)
end


# Generate Ntot_log_1d, Ftot_log_1d 

function generate_Ntot_grid(tab_r, tab_v, tab_m, rmin, rmax, vmin, vmax, dr, dv)
    Ni = Int(floor((log(rmax) - log(rmin)) / dr))
    Nj = Int(floor((log(vmax) - log(vmin)) / dv))
    m_sum_grid = zeros(Float64, Ni+1, Nj+1)
    
    for (r, v, m) in zip(tab_r, tab_v, tab_m)
        logr = log(r)
        logv = log(v)
        i = Int(floor((logr - log(rmin)) / dr))
        j = Int(floor((logv - log(vmin)) / dv))
        if 0 <= i <= Ni && 0 <= j <= Nj
            m_sum_grid[i+1, j+1] += m
        end
    end
    
    Ntot_log_grid = m_sum_grid ./ (dr * dv)
    return Ntot_log_grid, Ni, Nj
end

function Ntot_interpolator(Ntot_log_grid, Ni, Nj, rmin, rmax, vmin, vmax, dr, dv)
    function Ntot_log_1d(r, v)

        #Check bounds 
        if r < rmin || r > rmax || v < vmin || v > vmax
            return 0.0
        end
        
        logr = log(r)
        logv = log(v)

        # Bilinear interpolation
        il = clamp(Int(floor((logr - log(rmin)) / dr)), 0, Ni-1)
        jl = clamp(Int(floor((logv - log(vmin)) / dv)), 0, Nj-1)

        ir = il + 1
        jr = jl + 1
        
        ntot_ll = Ntot_log_grid[il+1, jl+1]
        ntot_lr = Ntot_log_grid[il+1, jr+1]
        ntot_rl = Ntot_log_grid[ir+1, jl+1]
        ntot_rr = Ntot_log_grid[ir+1, jr+1]
        
        x1 = log(rmin) + il * dr
        x2 = log(rmin) + ir * dr
        y1 = log(vmin) + jl * dv
        y2 = log(vmin) + jr * dv
        
        w11 = (x2 - logr) * (y2 - logv) / ((x2 - x1) * (y2 - y1))
        w12 = (x2 - logr) * (logv - y1) / ((x2 - x1) * (y2 - y1))
        w21 = (logr - x1) * (y2 - logv) / ((x2 - x1) * (y2 - y1))
        w22 = (logr - x1) * (logv - y1) / ((x2 - x1) * (y2 - y1))
        
        return w11 * ntot_ll + w12 * ntot_lr + w21 * ntot_rl + w22 * ntot_rr
    end
    
    return Ntot_log_1d
end

# Wrapper for grid, Ntot_log
function find_Ntot(tab_r, tab_v, tab_m, rmin, rmax, vmin, vmax, dr, dv)
    Ntot_log_grid, Ni, Nj = generate_Ntot_grid(tab_r, tab_v, tab_m, rmin, rmax, vmin, vmax, dr, dv)
    Ntot_log_1d = Ntot_interpolator(Ntot_log_grid, Ni, Nj, rmin, rmax, vmin, vmax, dr, dv)
    return Ntot_log_grid, Ntot_log_1d
end

function generate_Ftot_grid(tab_r, tab_v, tab_m, rmin, rmax, vmin, vmax, dr, dv)
    Ni = Int(floor((log(rmax) - log(rmin)) / dr))
    Nj = Int(floor((log(vmax) - log(vmin)) / dv))
    
    m2_sum_grid = zeros(Float64, Ni+1, Nj+1)
    
    for (r, v, m) in zip(tab_r, tab_v, tab_m)
        logr = log(r)
        logv = log(v)
        i = Int(floor((logr - log(rmin)) / dr))
        j = Int(floor((logv - log(vmin)) / dv))
        if 0 <= i <= Ni && 0 <= j <= Nj
            m2_sum_grid[i+1, j+1] += m ^ 2
        end
    end
    
    Ftot_log_grid = m2_sum_grid ./ (dr * dv)
    return Ftot_log_grid, Ni, Nj
end

function Ftot_interpolator(Ftot_log_grid, Ni, Nj, rmin, rmax, vmin, vmax, dr, dv)
    function Ftot_log_1d(r, v)
        # Check bounds
        if r < rmin || r > rmax || v < vmin || v > vmax
            return 0.0
        end
        
        logr = log(r)
        logv = log(v)

        # Bilinear interpolation
        il = clamp(Int(floor((logr - log(rmin)) / dr)), 0, Ni-1)
        jl = clamp(Int(floor((logv - log(vmin)) / dv)), 0, Nj-1)

        ir = il + 1
        jr = jl + 1
        
        ftot_ll = Ftot_log_grid[il+1, jl+1]
        ftot_lr = Ftot_log_grid[il+1, jr+1]
        ftot_rl = Ftot_log_grid[ir+1, jl+1]
        ftot_rr = Ftot_log_grid[ir+1, jr+1]
        
        x1 = log(rmin) + il * dr
        x2 = log(rmin) + ir * dr
        y1 = log(vmin) + jl * dv
        y2 = log(vmin) + jr * dv
        
        w11 = (x2 - logr) * (y2 - logv) / ((x2 - x1) * (y2 - y1))
        w12 = (x2 - logr) * (logv - y1) / ((x2 - x1) * (y2 - y1))
        w21 = (logr - x1) * (y2 - logv) / ((x2 - x1) * (y2 - y1))
        w22 = (logr - x1) * (logv - y1) / ((x2 - x1) * (y2 - y1))
        
        return w11 * ftot_ll + w12 * ftot_lr + w21 * ftot_rl + w22 * ftot_rr
    end
    
    return Ftot_log_1d
end

function find_Ftot(tab_r, tab_v, tab_m, rmin, rmax, vmin, vmax, dr, dv)
    Ftot_log_grid, Ni, Nj = generate_Ftot_grid(tab_r, tab_v, tab_m, rmin, rmax, vmin, vmax, dr, dv)
    Ftot_log_1d = Ftot_interpolator(Ftot_log_grid, Ni, Nj, rmin, rmax, vmin, vmax, dr, dv)
    return Ftot_log_grid, Ftot_log_1d
end

function read_in(filename)

    # read in data
    data=CSV.File(filename)
    omega_cen=DataFrame(data)

    # set up arrays
    #global id_tab, m_tab, r_tab, vr_tab, vt_tab, startype_tab, binflag = eachcol(omega_cen[:,1:7])
    global r_tab,m_tab,psi_tab,vr_tab,vt_tab=eachcol(omega_cen[:,1:5])
    global v_tab = sqrt.(vr_tab .^ 2 + vt_tab .^2)

    # find bounds
    global vmin=minimum(v_tab)
    global vmax=maximum(v_tab)

    global rmin=minimum(r_tab)
    global rmax=maximum(r_tab)

    # coulomb logarithm, test mass
    global coulomb_log = M_bh / mean(m_tab)
    global m_test = mean(m_tab) # Test mass

    # find potential  
    global U_tab, M_tab = find_psi_arrays(r_tab,m_tab)
    global psi_itp = psi_exact(U_tab,M_tab,r_tab)

    # calculate Ntot_log, Ftot_log
    global Ntot_grid,Ntot_log_1d=find_Ntot(r_tab,v_tab,m_tab,rmin,rmax,vmin,vmax,grid_size,grid_size)
    global Ftot_grid,Ftot_log_1d=find_Ftot(r_tab,v_tab,m_tab,rmin,rmax,vmin,vmax,grid_size,grid_size)
end

end