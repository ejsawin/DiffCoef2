module MakeGlobals

# --- Import necessary packages ---
using CSV
using DataFrames
using Interpolations
using Statistics

# --- Export variables, functions ---
export G, M_bh, N, m_test, coulomb_log
export dr, dv, int_steps

export m_tab, r_tab, vr_tab, vt_tab, v_tab
export rmin, rmax, vmin, vmax

export find_shell, find_psi_arrays, find_psi_exact, psi_exact
export psi_arr, M_arr, psi_calc

export get_sampling_DF_rv, bilinear_interp
export tabr_samp, tabv_samp, tabN_rv, tabN_rv_grid, tabF_rv, tabF_rv_grid
export N_1D, F_1D

export read_file 


# --- Physical constants ---
const G = 1
const M_bh = 0 #8.500250825e-03

# --- Numerical constants ---

int_steps=1000
dr=0.01
dv=0.01


# Potential Calculation

# --- Shell finding via bisection---
function find_shell(rr,r_tab)

    #check edge cases
    if rr < r_tab[1]
        return 0
    elseif rr > r_tab[N]
        return N
    elseif rr == r_tab[N]
        return N-1
    end

    #simple bisection
    left, right = 1, N - 1
    while left < right
        mid = div(left + right + 1, 2)
        if r_tab[mid] <= rr
            left = mid
        else
            right = mid - 1
        end
    end

    return left
end

# --- Generate psi, mass arrays via Henon with central BH ---
function find_psi_arrays(r_tab,m_tab)

    #set up array
    psi_arr=zeros(N+1)
    M_arr=zeros(N+1)

    psi_arr[N+1] = 0.0 # Eq 1.7
    M_arr[N] = sum(m_tab) # Eq 1.8
    psi_arr[N] = -G*(M_arr[N] + M_bh) / r_tab[N] #total enclosed mass

    #iterate over shells
    for k in (N-1):-1:1
        M_arr[k] = M_arr[k+1]-m_tab[k+1]
        psi_arr[k] = psi_arr[k+1]-G*M_arr[k]*(1/r_tab[k] - 1/r_tab[k+1]) - G * M_bh / r_tab[k] # Eq 1.15

    end

    #return arrays
    return psi_arr, M_arr
end

# --- Find exact potential using arrays ---
function find_psi_exact(rr,r_tab,psi_arr,M_arr)
    #find shell
    k=find_shell(rr,r_tab)

    #test if first, last shell
    if k == 0
        return psi_arr[1] - G * M_bh * (1/rr - 1/r_tab[1])
    
    elseif k == N
        return -G * (M_arr[N] + M_bh)/rr
    
    #between
    else 
        rk, rk1 = r_tab[k], r_tab[k+1]
        psi_k, psi_k1 = psi_arr[k], psi_arr[k+1]
        return psi_k+((1/rk - 1/rr)/(1/rk - 1/rk1))*(psi_k1-psi_k) # Eq 1.11
    end
end

# --- Wrapper function for psi ---
function psi_exact(r_tab,psi_arr,M_arr)
    return rr -> find_psi_exact(rr,r_tab,psi_arr,M_arr)
end

# Distribution function computation

# --- Find arrays, mass weighted distribution function ---
function get_sampling_DF_rv(tabr, tabv, tabm)

    # set up arrays from 0 -> (rmax-dr) 
    tabr_samp = collect(0.0:dr:rmax-dr) .+ dr/2 # bin edges -> center
    nbr = length(tabr_samp)
    tabv_samp = collect(0.0:dv:vmax-dv) .+ dv/2
    nbv = length(tabv_samp)

    tabF_rv = zeros(Float64, nbr, nbv) # grid
    tabF_rv_grid = zeros(Float64, nbr * nbv) # array

    tabN_rv = zeros(Float64, nbr, nbv) # grid
    tabN_rv_grid = zeros(Float64, nbr * nbv) # array

    Threads.@threads for i=1:N 

        # r = (ir-1)*dr + deltar
        # v = (iv-1)*dv + deltav
        # index-1 = iv-1 + (ir-1)*nbv

        r = tabr[i]
        v = tabv[i]
        mass = tabm[i]
        
        # bin indices
        ir = floor(Int64, r/dr) + 1
        iv = floor(Int64, v/dv) + 1

        # indices for array 
        index = iv + (ir-1)*nbv

        # check within bin, sum masses
        if ((ir <= nbr) && (iv <= nbv))

            tabN_rv[ir, iv] += mass/(dr*dv)
            tabN_rv_grid[index] += mass/(dr*dv)

            tabF_rv[ir, iv] += (mass^2)/(dr*dv)
            tabF_rv_grid[index] += (mass^2)/(dr*dv)
        end

    end

    return tabr_samp, tabv_samp, tabN_rv, tabN_rv_grid, tabF_rv, tabF_rv_grid
end


# --- Bilinear interpolation (weighted means) for grid -> integrals ---
function bilinear_interp(rr, vv, rgrid, vgrid, fgrid)
    # assuming rgrid, vgrid are centers
    dr = rgrid[2] - rgrid[1]
    dv = vgrid[2] - vgrid[1]

    # find grid edges
    r_edges = [rgrid[1]-dr/2; rgrid .+ dr/2]
    v_edges = [vgrid[1]-dv/2; vgrid .+ dv/2]

    # find lower-left cell
    ir = searchsortedlast(r_edges, rr)
    iv = searchsortedlast(v_edges, vv)

    # check bounds
    if ir < 1 || iv < 1 || ir >= size(fgrid,1) || iv >= size(fgrid,2)
        return 0.0
    end

    # coordinates of corners
    r1, r2 = r_edges[ir], r_edges[ir+1]
    v1, v2 = v_edges[iv], v_edges[iv+1]

    # function values
    f11 = fgrid[ir, iv]
    f12 = fgrid[ir, iv+1]
    f21 = fgrid[ir+1, iv]
    f22 = fgrid[ir+1, iv+1]

    denom = (r2 - r1) * (v2 - v1)

    # weights
    w11 = (r2 - rr) * (v2 - vv) / denom
    w12 = (r2 - rr) * (vv - v1) / denom
    w21 = (rr - r1) * (v2 - vv) / denom
    w22 = (rr - r1) * (vv - v1) / denom

    return w11*f11 + w12*f12 + w21*f21 + w22*f22
end

# --- Read in data, compute psi and DFs ---
function read_file(filename)

    #read in data, make DataFrame
    data=CSV.File(filename)
    data_frame=DataFrame(data)

    #set up arrays, calculate velocity
    #global id_tab, m_tab, r_tab, vr_tab, vt_tab, startype_tab, binflag = eachcol(omega_cen[:,1:7])
    global r_tab, m_tab, psi_tab, vr_tab, vt_tab = eachcol(data_frame[:,1:5])
    global v_tab = sqrt.(vr_tab .^ 2 + vt_tab .^2)

    # Determine array length
    global N = length(r_tab)

    #find (r,v) bounds
    global rmin, rmax = extrema(r_tab)
    global vmin, vmax = extrema(v_tab)

    #calculate coulomb logarithm, test mass
    global coulomb_log = log(0.11*N)# Uniform mass- usually M_bh / mean(m_tab)
    global m_test = mean(m_tab)

    #find psi, m arrays -> calculate potential
    global psi_arr, M_arr = find_psi_arrays(r_tab,m_tab)
    global psi_calc = psi_exact(r_tab,psi_arr,M_arr)

    #compute DataFrame arrays
    global tabr_samp, tabv_samp, tabN_rv, tabN_rv_grid, tabF_rv, tabF_rv_grid = get_sampling_DF_rv(r_tab,v_tab,m_tab)

    #compute distribution function
    global N_1D = (r,v) -> bilinear_interp(r,v,tabr_samp,tabv_samp,tabN_rv)
    global F_1D = (r,v) -> bilinear_interp(r,v,tabr_samp,tabv_samp,tabF_rv)

end

end