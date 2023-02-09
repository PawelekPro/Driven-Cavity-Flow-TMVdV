function [output_u, output_v] = boundary_conditions(Lx, Ly, Nx, Ny, u, v, time)
f_time = time;
bc_top_wall = linspace(0,Lx,Nx);
bt = f_time*(1 - cos(2*pi*bc_top_wall));

u(Nx,:) = 0;
u(:,1) = 0;
u(:,Ny) = 0;
u(1:Ny,1) = bt(1:Nx);

v(Nx,:) = 0;
v(:,1) = 0;
v(:,Ny) = 0;
v(1,:) = 0;

output_u = u;
output_v = v;
end

