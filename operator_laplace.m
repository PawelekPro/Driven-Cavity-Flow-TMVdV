function output = operator_laplace(Nx, Ny, dx, dy)
%Lux
ud = diag(ones(1,Nx*Ny-1),1);
bd = diag(ones(1,Nx*Ny-1),-1);
Lux = (-2 * eye(Nx*Ny) + ud + bd);

%Luy
ud = diag(ones(1,Nx*Nx-Ny),Ny);
bd = diag(ones(1,Nx*Nx-Ny),-Ny);
Luy = (-2 * eye(Nx*Ny) + ud + bd);

output = 1/(dx^2) * Lux + 1/(dy^2) * Luy;




end