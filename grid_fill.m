function grid_nodes = grid_fill(Nx, Ny, Lx, Ly)
% 1 = interior
% 2 = boundary top, 
% 3 = boundary right, 
% 4 = boundary left, 
% 5 = boundary bot

grid_nodes = struct('number',{},...
    'type',{},...
    'code',{},...
    'coord_x',{},...
    'coord_y',{});

for i = 1:Nx*Ny
    grid_nodes(i).code = 1;
    grid_nodes(i).type = 'int';
    grid_nodes(i).number = i;
end

for i = 1 : Nx : Nx*Ny
    grid_nodes(i).type = 'bl';
    grid_nodes(i).code = 4;
end

for i = Nx : Nx : Nx*Ny
    grid_nodes(i).type = 'br';
    grid_nodes(i).code = 3;
end

for i = 1:Nx
    grid_nodes(i).type = 'bt';
    grid_nodes(i).code = 2;
end

for i = (Nx*(Ny-1)+1):Nx*Ny
    grid_nodes(i).type = 'bb';
    grid_nodes(i).code = 5;
end


xce = linspace(0,Lx,Nx);
yce = linspace(0,Ly,Ny);
[Xce,Yce] = meshgrid(xce,yce);
k = 1;
for i = 1:(Ny)
    for j = 1:(Nx)
        grid_nodes(k).coord_x = Xce(i,j);
        grid_nodes(k).coord_y = Xce(j,i);
        k = k + 1;
    end
end
end


