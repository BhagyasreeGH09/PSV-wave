using LinearAlgebra, Plots

Vp = 4000.0              
nu = 0.25               
rho = 2500.0            
Vs = Vp * sqrt((0.5 - nu) / (1 - nu))  
λ = rho * (Vp^2 - 2*Vs^2)             
μ = rho * Vs^2                        

Nx, Nz = 601,601          
dx, dz = 25,25    
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

p1, p2, p3, p4 = Int(1500/dx), Int(2000/dx), Int(2500/dx), Int(3000/dx)
ir1, ir2, ir3, ir4 = isrc + p1, isrc + p2, isrc + p3, isrc + p4
jr = 2

u_r1 = zeros(nt)
u_r2 = zeros(nt)
u_r3 = zeros(nt)
u_r4 = zeros(nt)

for n in 1:nt

    Txx[isrc, jsrc] += source_time[n]
    Tzz[isrc, jsrc] += source_time[n]

    for i in 2:Nx-1, j in 2:Nz-1
        dTxx_dx = (Txx[i,j] - Txx[i-1,j]) / dx
        dTxz_dz = (Txz[i,j] - Txz[i,j-1]) / dz
        Vx[i,j] += dt / rho * (dTxx_dx + dTxz_dz)
    end

    for i in 1:Nx-1, j in 1:Nz-1
        dTxz_dx = (Txz[i+1,j] - Txz[i,j]) / dx
        dTzz_dz = (Tzz[i,j+1] - Tzz[i,j]) / dz
        Vz[i,j] += dt / rho * (dTxz_dx + dTzz_dz)
    end

    for j in 1:Nz
        Vx[1,j] = Vx[1,j] + ((dt*Vp)/dx) * (Vx[2,j] - Vx[1,j])
        Vx[end,j] = Vx[end,j] + ((dt*Vp)/dx)* (Vx[end,j] - Vx[end-1,j])
    end

    for i in 2:Nx-1
        Vx[i,end] = Vx[i,end] + ((dt*Vp)/dz)* (Vx[i,end] - Vx[i,end-1])
    end

    for i in 1:Nx-1, j in 2:Nz-1
        dVx_dx = (Vx[i+1,j] - Vx[i,j]) / dx
        dVz_dz = (Vz[i,j] - Vz[i,j-1]) / dz
        Txx[i,j] += dt * ((λ + 2μ) * dVx_dx + λ * dVz_dz)
        Tzz[i,j] += dt * (λ * dVx_dx + (λ + 2μ) * dVz_dz)
    end

    for i in 2:Nx-1, j in 1:Nz-1
       dVz_dx = (Vz[i,j] - Vz[i-1,j]) / dx
       dVx_dz = (Vx[i,j+1] - Vx[i,j]) / dz
       Txz[i,j] += dt * μ * (dVz_dx + dVx_dz)
    end

    for i in 1:Nx-1
        Txx[i,end] = Txx[i,end] + ((dt*Vp)/dz)* (Txx[i,end] - Txx[i,end-1])
        Tzz[i,end] = Tzz[i,end] + ((dt*Vp)/dz)* (Tzz[i,end] - Tzz[i,end-1])
    end
    for j in 1:Nz-1
        Txz[1,j] = Txz[1,j] + ((dt*Vp)/dx) * (Txz[2,j] - Txz[1,j])
        Txz[end,j] = Txz[end,j] + ((dt*Vp)/dx)* (Txz[end,j] - Txz[end-1,j])
    end

    Tzz[:,1] .= 0.0
    Txz[:,1] .= -Txz[:,2] 

 u_r1[n] = Vx[ir1, jr]
 u_r2[n] = Vx[ir2, jr]
 u_r3[n] = Vx[ir3, jr]
 u_r4[n] = Vx[ir4, jr]
   
end
u_disp1 = cumsum(u_r1) .* dt
u_disp2 = cumsum(u_r2) .* dt
u_disp3 = cumsum(u_r3) .* dt
u_disp4 = cumsum(u_r4) .* dt

u_disp1 = u_disp1 / maximum(abs.(u_disp1))
u_disp2 = u_disp2 / maximum(abs.(u_disp2))
u_disp3 = u_disp3 / maximum(abs.(u_disp3))
u_disp4 = u_disp4 / maximum(abs.(u_disp4))

plot(t, u_disp1, lw=2, xlabel="Time (s)", ylabel="Normalized radial displacement",
    legend=false)
savefig("seismogram_1500m.png")

plot(t, u_disp2, lw=2, xlabel="Time (s)", ylabel="Normalized radial displacement",
    legend=false)
savefig("seismogram_2000m.png")

plot(t, u_disp3, lw=2, xlabel="Time (s)", ylabel="Normalized radial displacement",
      legend=false)
savefig("seismogram_2500m.png")

plot(t, u_disp4, lw=2, xlabel="Time (s)", ylabel="Normalized radial displacement",
     legend=false)
savefig("seismogram_3000m.png")