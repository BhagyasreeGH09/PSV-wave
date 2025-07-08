using LinearAlgebra, Plots

Vp = 4000.0              
nu = 0.25               
rho = 2500.0            
Vs = Vp * sqrt((0.5 - nu) / (1 - nu))  
λ = rho * (Vp^2 - 2*Vs^2)             
μ = rho * Vs^2                        

Nx, Nz = 201,201          
dx, dz = 50,50    
dt = 0.4 * dx / (Vp * sqrt(2))
T = 2
nt = Int(round(T/dt))     
t = collect(0:dt:(nt-1)*dt)
           
a = 200
t0=0.1
#source_time = -2 .* a .* (t .- t0) .* exp.(-a .* (t .- t0).^2)
source_time =exp.(-a .* (t .- t0).^2)
Txx = zeros(Nx-1, Nz)
Tzz = zeros(Nx-1, Nz)
Txz = zeros(Nx, Nz-1)
Vx = zeros(Nx, Nz)
Vz = zeros(Nx-1, Nz-1)

isrc, jsrc = div(Nx-1,2), 1

p=40
ir, jr = isrc + p, jsrc
u_r = zeros(Float64, nt)

function apply_clayton_enquist!(Vx, Vz, Vp, dt, dx, dz)
    for j in 2:Nz-2
        Vx[1,j] = Vx[1,j] + ((dt*Vp)/dx) * (Vx[2,j] - Vx[1,j])
        Vx[end,j] = Vx[end,j] - ((dt*Vp)/dx)* (Vx[end,j] - Vx[end-1,j])
    end

    for i in 2:Nx-2
        Vx[i,end] = Vx[i,end] - ((dt*Vp)/dz)* (Vx[i,end] - Vx[i,end-1])
    end
end

for n in 1:nt

    Txx[isrc, jsrc] += source_time[n]
    Tzz[isrc, jsrc] += source_time[n]

    for i in 2:Nx-1, j in 1:Nz-1
        dTxx_dx = (Txx[i,j] - Txx[i-1,j]) / dx
        dTxz_dz = (Txz[i,j+1] - Txz[i,j]) / dz
        Vx[i,j] += dt / rho * (dTxx_dx + dTxz_dz)
    end

    for i in 2:Nx-2, j in 1:Nz-2
        dTxz_dx = (Txz[i+1,j] - Txz[i,j]) / dx
        dTzz_dz = (Tzz[i,j+1] - Tzz[i,j]) / dz
        Vz[i,j] += dt / rho * (dTxz_dx + dTzz_dz)
    end

   apply_clayton_enquist!(Vx, Vz, Vp, dt, dx, dz)


    for i in 1:Nx-1, j in 1:Nz-1
        dVx_dx = (Vx[i+1,j] - Vx[i,j]) / dx
        dVz_dz = (Vz[i,j+1] - Vz[i,j]) / dz
        Txx[i,j] += dt * ((λ + 2μ) * dVx_dx + λ * dVz_dz)
        Tzz[i,j] += dt * (λ * dVx_dx + (λ + 2μ) * dVz_dz)
    end

    for i in 2:Nx, j in 2:Nz-1
     dVz_dx = (Vz[i,j] - Vz[i-1,j]) / dx
     dVx_dz = (Vx[i,j] - Vx[i,j-1]) / dz
     Txz[i,j] += dt * μ * (dVz_dx + dVx_dz)
    end

    Tzz[:,1] .= 0.0
    Txz[:,1] .= -Txz[:,2] 

    u_r[n] = Vx[ir, jr]
   
end
u_disp = cumsum(u_r) .* dt

plot(t, u_disp ./ maximum(abs.(u_disp)), lw=2, xlabel="Time (s)", ylabel="Normalized displacement", title="2D Elastic Wave - Displacement Seismogram at x=$(p*dx)", legend=false)
# savefig("Gaussian_freesurface.png")