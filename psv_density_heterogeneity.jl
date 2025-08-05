using LinearAlgebra, Plots

Vp = 4000.0   
Nx, Nz = 201, 201          
dx, dz = 50, 50       
dt = 0.4 * dx / (Vp * sqrt(2))   
T = 2
nt = Int(round(T/dt))     
t = collect(0:dt:(nt-1)*dt)
           
nu = 0.25  
Vs = Vp * sqrt((0.5 - nu) / (1 - nu))             
rho= fill(3500.0, Nx, Nz) 
λ = zeros(Nx,Nz)
μ = zeros(Nx,Nz)
for i in 1:Nx, j in 1:Nz
    if j < Nz ÷ 4
        rho[i,j] = 2000.0
    end
   λ[i,j] = rho[i,j]* (Vp^2 - 2*Vs^2)  
   μ[i,j] = rho[i,j] * Vs^2 
end         
          
a = 200
t0=0.1
source_time=-2*a.*(t.-t0).*exp.(-a .* (t .- t0).^2)

Txx = zeros(Nx-1, Nz)
Tzz = zeros(Nx-1, Nz)
Txz = zeros(Nx, Nz-1)
Vx = zeros(Nx, Nz)
Vz = zeros(Nx-1, Nz-1)

isrc, jsrc = div(Nx-1,2), div(Nz-1,2)

times_to_plot = [0.5, 1.0, 1.4,1.6,1.8,2]  
frames_to_plot = [Int(round(ti / dt)) for ti in times_to_plot]

x = collect(0:dx:(Nx-2)*dx) 
z = collect(0:dz:(Nz-1)*dz)

for n in 1:nt

    Txx[isrc, jsrc] += source_time[n]
    Tzz[isrc, jsrc] += source_time[n]

   for i in 2:Nx-1, j in 2:Nz-1
        dTxx_dx = (Txx[i,j] - Txx[i-1,j]) / dx
        dTxz_dz = (Txz[i,j] - Txz[i,j-1]) / dz
        Vx[i,j] += dt / rho[i,j] * (dTxx_dx + dTxz_dz)
   end

    for i in 1:Nx-2, j in 1:Nz-1
        dTxz_dx = (Txz[i+1,j] - Txz[i,j]) / dx
        dTzz_dz = (Tzz[i,j+1] - Tzz[i,j]) / dz
        Vz[i,j] += dt / rho[i,j] * (dTxz_dx + dTzz_dz)
    end

    for j in 1:Nz
        Vx[1,j] = Vx[1,j] + ((dt*Vp)/dx) * (Vx[2,j] - Vx[1,j]) 
        Vx[end,j] = Vx[end,j] - ((dt*Vp)/dx)* (Vx[end,j] - Vx[end-1,j])
    end

    for i in 2:Nx-1
        Vx[i,end] = Vx[i,end] - ((dt*Vs)/dz)* (Vx[i,end] - Vx[i,end-1])
        Vx[i,1] = Vx[i,1] +((dt*Vs)/dz)* (Vx[i,2] - Vx[i,1])
    end

    for i in 1:Nx-1, j in 2:Nz-1
        dVx_dx = (Vx[i+1,j] - Vx[i,j]) / dx
        dVz_dz = (Vz[i,j] - Vz[i,j-1]) / dz
        Txx[i,j] += dt * ((λ[i,j] + 2μ[i,j]) * dVx_dx + λ[i,j] * dVz_dz)
        Tzz[i,j] += dt * (λ[i,j] * dVx_dx + (λ[i,j] + 2μ[i,j]) * dVz_dz)
    end

    for i in 2:Nx-1, j in 1:Nz-1
       dVz_dx = (Vz[i,j] - Vz[i-1,j]) / dx
       dVx_dz = (Vx[i,j+1] - Vx[i,j]) / dz
       Txz[i,j] += dt * μ[i,j] * (dVz_dx + dVx_dz)
    end 

 
    if n in frames_to_plot
        heatmap(x,z,Txx', color=:seismic, clims=(-1e-3, 1e-3), xlabel="Range in meters", ylabel="Range in meters",
             framestyle=:box)
        savefig("Derivative frame3_t$(round(n*dt, digits=3)).png")  
    end
   
end
