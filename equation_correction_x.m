function [LHS, RHS] = equation_correction_x(Lx, Ly, Nx, Ny, u, v, nodes, LHS, RHS, time)

[u, v] = boundary_conditions(Lx,Ly,Nx,Ny,u,v,time);

for i = 1:Nx*Ny
    if(nodes(i).code == 2) % top wall
        LHS(i, nodes(i).number) = 1;
        RHS(nodes(i).number,1) = u(nodes(i).number,1); 
    elseif(nodes(i).code == 4) % left wall
        LHS(i, nodes(i).number) = 1;
        RHS(nodes(i).number,1) = 0;
    elseif(nodes(i).code == 3) % right wall
        LHS(i, nodes(i).number) = 1;
        RHS(nodes(i).number,1) = 0;
    elseif(nodes(i).code == 5) % bottom wall
        LHS(i, nodes(i).number) = 1;
        RHS(nodes(i).number,1) = 0;
    end
end

end