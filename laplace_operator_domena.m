function output = laplace_operator_domena(Nx, Ny, dx, dy, nodes)

% kierunek x
top_boundary_nodes_x = zeros(Nx*Ny);
left_boundary_nodes_X = zeros(Nx*Ny);
right_boundary_nodes_x = zeros(Nx*Ny);
bot_boundary_nodes_x = zeros(Nx*Ny);
interior_nodes_x = zeros(Nx*Ny);

for i = 1:Nx*Ny
    if(nodes(i).code == 2)
        top_boundary_nodes_x(i,nodes(i).number) = 1;
    elseif(nodes(i).code == 4)
        left_boundary_nodes_X(i,nodes(i).number) = 1;
    elseif(nodes(i).code == 3)
        right_boundary_nodes_x(i,nodes(i).number) = 1;
    elseif(nodes(i).code == 5)
        bot_boundary_nodes_x(i,nodes(i).number) = 1;
    elseif(nodes(i).code == 1)
        interior_nodes_x(i,nodes(i).number) = -2;
        interior_nodes_x(i,nodes(i).number - 1) = 1;
        interior_nodes_x(i,nodes(i).number + 1) = 1;
    end
end


% kierunek y
interior_nodes_y = zeros(Nx*Ny);

for i = 1:Nx*Ny
    if(nodes(i).code == 1)
        interior_nodes_y(i,nodes(i).number) = -2;
        interior_nodes_y(i,nodes(i).number - Nx) = 1;
        interior_nodes_y(i,nodes(i).number + Nx) = 1;
    end
end

output = (1/dy^2)*interior_nodes_y + (1/dx^2)*interior_nodes_x + bot_boundary_nodes_x + top_boundary_nodes_x + left_boundary_nodes_X + right_boundary_nodes_x;
output = sparse(output);
end
