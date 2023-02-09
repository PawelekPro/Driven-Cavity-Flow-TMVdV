function [px,py] = grad_pressure(Nx, Ny, dx, dy)

% CENTRAL DIFFERENCES
u_xd = diag(ones(1,Nx*Ny-1),1);
u_xd(Nx:Nx:end) = 0;

b_xd = -diag(ones(1,Nx*Ny-1),-1);
b_xd((Nx+1):Nx:end) = 0;

px = 1/(2*dx)*(u_xd + b_xd);


u_yd = diag(ones(1,Nx*Ny-Ny),Nx);
b_yd = -diag(ones(1,Nx*Ny-Ny),-Nx);

py = 1/(2*dy)*(u_yd + b_yd);


sparse(px);
sparse(py);
end