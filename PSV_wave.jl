using LinearAlgebra, Plots

Vp = 4000.0              
nu = 0.25               
rho = 2500.0            
Vs = Vp * sqrt((0.5 - nu) / (1 - nu))  
λ = rho * (Vp^2 - 2*Vs^2)             
μ = rho * Vs^2                        

Nx, Nz = 201, 201          
dx, dz = 50, 50       
dt = 0.4 * dx / (Vp * sqrt(2))   
T = 2
nt = Int(round(T/dt))     
t = collect(0:dt:(nt-1)*dt)
           
a = 200
t0=0.1
#source_time =exp.(-a .* (t .- t0).^2)
#source_time = (1 .- 2 .* (π .* a .* (t .- t0)).^2) .* exp.(-(π .* a .* (t .- t0)).^2)
source_time=-2*a.*(t.-t0).*exp.(-a .* (t .- t0).^2)
Txx = zeros(Nx-1, Nz-1)
Tzz = zeros(Nx-1, Nz-1)
Txz = zeros(Nx, Nz)
Vx = zeros(Nx, Nz-1)
Vz = zeros(Nx-1, Nz)

isrc, jsrc = div(Nx-1,2), div(Nz-1,2)

ir, jr = isrc + 20, jsrc
u_r = zeros(Float64, nt)

times_to_plot = [0.5, 1.0, 1.4,1.6,1.8,2]   # seconds
frames_to_plot = [Int(round(ti / dt)) for ti in times_to_plot]

function apply_clayton_enquist!(Vx, Vz, Vp, dt, dx, dz)
 for j in 2:Nz-2
        Vx[1,j] = Vx[1,j] + ((dt*Vp)/dx) * (Vx[2,j] - Vx[1,j])
        Vz[1,j] = Vz[1,j] + ((dt*Vs)/dx) * (Vz[2,j] - Vz[1,j])
    end  

    for j in 2:Nz-2
        Vx[end,j] = Vx[end,j] - ((dt*Vp)/dx)* (Vx[end,j] - Vx[end-1,j])
        Vz[end,j] = Vz[end,j] - ((dt*Vs)/dx)* (Vz[end,j] - Vz[end-1,j])  
    end

    for i in 2:Nx-2
        Vx[i,end] = Vx[i,end] - ((dt*Vs)/dz)* (Vx[i,end] - Vx[i,end-1])
        Vz[i,end] = Vz[i,end] - ((dt*Vp)/dz) * (Vz[i,end] - Vz[i,end-1])
    end

     for i in 2:Nx-2
        Vx[i,1] = Vx[i,1] +((dt*Vs)/dz)* (Vx[2,1] - Vx[1,i])
        Vz[i,1] = Vz[i,1] + ((dt*Vp)/dz) * (Vz[2,i] - Vz[1,i])
    end
end
for n in 1:nt
   for i in 2:Nx-1, j in 1:Nz-2
        dTxx_dx = (Txx[i,j] - Txx[i-1,j]) / dx
        dTxz_dz = (Txz[i,j+1] - Txz[i,j]) / dz
        Vx[i,j] += dt / rho * (dTxx_dx + dTxz_dz)
   end

    for i in 2:Nx-2, j in 2:Nz-1
        dTxz_dx = (Txz[i+1,j] - Txz[i,j]) / dx
        dTzz_dz = (Tzz[i,j] - Tzz[i,j-1]) / dz
        Vz[i,j] += dt / rho * (dTxz_dx + dTzz_dz)
    end

   #apply_clayton_enquist!(Vx, Vz, Vp, dt, dx, dz)

    for i in 2:Nx-2, j in 2:Nz-2
        dVx_dx = (Vx[i+1,j] - Vx[i,j]) / dx
        dVz_dz = (Vz[i,j+1] - Vz[i,j]) / dz
        Txx[i,j] += dt * ((λ + 2μ) * dVx_dx + λ * dVz_dz)
        Tzz[i,j] += dt * (λ * dVx_dx + (λ + 2μ) * dVz_dz)
    end

    for i in 2:Nx-1, j in 2:Nz-1
     dVz_dx = (Vz[i,j] - Vz[i-1,j]) / dx
     dVx_dz = (Vx[i,j] - Vx[i,j-1]) / dz
     Txz[i,j] += dt * μ * (dVz_dx + dVx_dz)
    end  

    Txx[isrc, jsrc] += source_time[n]
    Tzz[isrc, jsrc] += source_time[n]

    u_r[n] = Vx[ir, jr]
     if n in frames_to_plot
        heatmap(Txx', color=:seismic, clims=(-1e-3, 1e-3), xlabel="x", ylabel="z",
            title="Gaussian Derivative wavefield at t=$(round(n*dt, digits=3)) s", framestyle=:box)
        savefig("Derivative frame_t$(round(n*dt, digits=3)).png")  # save each plot if needed
    end
   
end
