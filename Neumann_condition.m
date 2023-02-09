function output = Neumann_condition(Nx,Ny,Lx,Ly,nodes)

dx = Lx/Nx;
dy = Ly/Ny;

% kierunek x
top_boundary_nodes_x = zeros(Nx*Ny);
left_boundary_nodes_x = zeros(Nx*Ny);
right_boundary_nodes_x = zeros(Nx*Ny);
bot_boundary_nodes_x = zeros(Nx*Ny);
interior_nodes_x = zeros(Nx*Ny);

for i = 1:Nx*Ny
    if(nodes(i).code == 2)
        if(nodes(i).number == 1)
            top_boundary_nodes_x(i,nodes(i).number) = -2;
            top_boundary_nodes_x(i,nodes(i).number + 1) = 1;%
        elseif(nodes(i).number == Nx)
            top_boundary_nodes_x(i,nodes(i).number) = -2;
            top_boundary_nodes_x(i,nodes(i).number - 1) = 1;%
        else
            top_boundary_nodes_x(i,nodes(i).number) = -2;
            top_boundary_nodes_x(i,nodes(i).number + 1) = 1;
            top_boundary_nodes_x(i,nodes(i).number - 1) = 1;
        end
    elseif(nodes(i).code == 4)
        left_boundary_nodes_x(i,nodes(i).number) = -2;
        left_boundary_nodes_x(i,nodes(i).number + 1) = 2;
        left_boundary_nodes_x(i,nodes(i).number - 1) = 0;
    elseif(nodes(i).code == 3)
        right_boundary_nodes_x(i,nodes(i).number) = -2;
        right_boundary_nodes_x(i,nodes(i).number + 1) = 0;
        right_boundary_nodes_x(i,nodes(i).number - 1) = 2;
    elseif(nodes(i).code == 5)
        if(nodes(i).number == Nx*Ny - Nx + 1)
            bot_boundary_nodes_x(i,nodes(i).number) = -2;
            bot_boundary_nodes_x(i,nodes(i).number - 1) = 0;
            bot_boundary_nodes_x(i,nodes(i).number + 1) = 2;
        elseif(nodes(i).number == Nx*Ny)
            bot_boundary_nodes_x(i,nodes(i).number) = -2;
            bot_boundary_nodes_x(i,nodes(i).number - 1) = 2;
        else
            bot_boundary_nodes_x(i,nodes(i).number) = -2;
            bot_boundary_nodes_x(i,nodes(i).number - 1) = 1;
            bot_boundary_nodes_x(i,nodes(i).number + 1) = 1;
        end   
    elseif(nodes(i).code == 1)
        interior_nodes_x(i,nodes(i).number) = -2;
        interior_nodes_x(i,nodes(i).number - 1) = 1;
        interior_nodes_x(i,nodes(i).number + 1) = 1;
    end
end

direct_x = 1/(dx^2) * (top_boundary_nodes_x + left_boundary_nodes_x + right_boundary_nodes_x + bot_boundary_nodes_x + interior_nodes_x);


% kierunek y
top_boundary_nodes_y = zeros(Nx*Ny);
left_boundary_nodes_y = zeros(Nx*Ny);
right_boundary_nodes_y = zeros(Nx*Ny);
bot_boundary_nodes_y = zeros(Nx*Ny);
interior_nodes_y = zeros(Nx*Ny);

for i = 1:Nx*Ny
    if(nodes(i).code == 1)
        interior_nodes_y(i,nodes(i).number) = -2;
        interior_nodes_y(i,nodes(i).number - Nx) = 1;
        interior_nodes_y(i,nodes(i).number + Nx) = 1;
    elseif(nodes(i).code == 2)
        if(nodes(i).number == 1 || nodes(i).number == Nx)  
            top_boundary_nodes_y(i,nodes(i).number) = -2;
            top_boundary_nodes_y(i,nodes(i).number + Nx) = 2;
        else
            top_boundary_nodes_y(i,nodes(i).number) = -2;
            top_boundary_nodes_y(i,nodes(i).number + Nx) = 2;
        end
    elseif(nodes(i).code == 4)
        left_boundary_nodes_y(i, nodes(i).number) = -2;
        left_boundary_nodes_y(i, nodes(i).number + Nx) = 1;
        left_boundary_nodes_y(i, nodes(i).number - Nx) = 1;
    elseif(nodes(i).code == 3)
        right_boundary_nodes_y(i, nodes(i).number) = -2;
        right_boundary_nodes_y(i, nodes(i).number + Nx) = 1;
        right_boundary_nodes_y(i, nodes(i).number - Nx) = 1;
    elseif(nodes(i).code == 5)
        bot_boundary_nodes_y(i, nodes(i).number) = -2;
        bot_boundary_nodes_y(i, nodes(i).number - Nx) = 2;
    end
end
direct_y = 1/(dy^2) * (top_boundary_nodes_y + left_boundary_nodes_y + right_boundary_nodes_y + bot_boundary_nodes_y + interior_nodes_y);

output = direct_x + direct_y;
output = sparse(output);
end


