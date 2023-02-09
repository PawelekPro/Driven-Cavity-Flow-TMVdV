function [uUx,vUy] = NS_term1(Nx, Ny, dx, dy, nodes, vel_u, vel_v)
uUx = zeros(Nx*Ny);
vUy = zeros(Nx*Ny);

u_xd = diag(ones(1,Nx*Nx-1),1);
u_xd(Nx:Nx:end) = 0;

b_xd = -diag(ones(1,Nx*Ny-1),-1);
b_xd((Nx+1):Nx:end) = 0;

d_x = 1/(2*dx)*(u_xd + b_xd);


u_yd = diag(ones(1,Nx*Ny-Ny),Nx);
b_yd = -diag(ones(1,Nx*Ny-Ny),-Nx);

d_y = 1/(2*dy)*(u_yd + b_yd);

for i = 1:Nx*Ny
   uUx(nodes(i).number,:) = vel_u(nodes(i).number,1) .* d_x(nodes(i).number,:);
   vUy(nodes(i).number,:) = vel_v(nodes(i).number,1) .* d_y(nodes(i).number,:);
end

uUx = sparse(uUx);
vUy = sparse(vUy);
end