% UPWIND
% function [output] = NS_term2(Nx,Ny,dx,dy,vel_u,vel_v,nodes)
% 
% output_x = zeros(Nx*Ny);
% output_y = zeros(Nx*Ny);
% 
%     for i = 1:Nx*Ny
%         if(nodes(i).code == 1)
%             output_x(nodes(i).number,i) = (vel_u(nodes(i).number,1) - vel_u(nodes(i).number - 1,1));
%             output_y(nodes(i).number,i) = (vel_v(nodes(i).number,1) - vel_v(nodes(i).number - Nx,1));
%         elseif(nodes(i).code == 2)
%             if(nodes(i).number == 1 )
%                 output_x(nodes(i).number,i) = (vel_u(nodes(i).number,1));
%                 output_y(nodes(i).number,i) = (vel_v(nodes(i).number,1));
%             elseif(nodes(i).number == Nx)
%                 output_x(nodes(i).number,i) = (vel_u(nodes(i).number,1) - vel_u(nodes(i).number-1,1));
%                 output_y(nodes(i).number,i) = (vel_v(nodes(i).number,1));
%             else
%                 output_x(nodes(i).number,i) = (vel_u(nodes(i).number,1) - vel_u(nodes(i).number - 1,1));
%                 output_y(nodes(i).number,i) = (vel_v(nodes(i).number,1));
%             end
%         
%         elseif(nodes(i).code == 5)
%             if(nodes(i).number == Nx*Ny - Nx + 1)
%                 output_x(nodes(i).number,i) = (vel_u(nodes(i).number,1));
%                 output_y(nodes(i).number,i) = (vel_v(nodes(i).number,1) - vel_v(nodes(i).number - Nx,1));
%             elseif(nodes(i).number == Nx*Ny)
%                 output_x(nodes(i).number,i) = (vel_u(nodes(i).number,1) - vel_u(nodes(i).number-1,1));
%                 output_y(nodes(i).number,i) = (vel_v(nodes(i).number,1) - vel_v(nodes(i).number - Nx,1));
%             else
%                 output_x(nodes(i).number,i) = (vel_u(nodes(i).number,1) - vel_u(nodes(i).number - 1,1));
%                 output_y(nodes(i).number,i) = (vel_v(nodes(i).number,1) - vel_v(nodes(i).number - Nx,1));
%             end
%         elseif(nodes(i).code == 4)
%             output_x(nodes(i).number,i) = (vel_u(nodes(i).number,1));
%             output_y(nodes(i).number,i) = (vel_v(nodes(i).number,1) - vel_v(nodes(i).number - Nx,1));
%         elseif(nodes(i).code == 3)
%             output_x(nodes(i).number,i) = (vel_u(nodes(i).number,1) - vel_u(nodes(i).number - 1,1));
%             output_y(nodes(i).number,i) = (vel_v(nodes(i).number,1) - vel_v(nodes(i).number - Nx,1));
%         end
%     end
%     output = 1/dx*output_x + 1/dy*output_y;
%     output = sparse(output);
% end

%CENTRAL DIFFERENCES
function [output] = NS_term2(Nx,Ny,dx,dy,vel_u,vel_v,nodes)

output_x = zeros(Nx*Ny);
output_y = zeros(Nx*Ny);

    for i = 1:Nx*Ny
        if(nodes(i).code == 1)
            output_x(nodes(i).number,i) = (vel_u(nodes(i).number + 1,1) - vel_u(nodes(i).number - 1,1));
            output_y(nodes(i).number,i) = (vel_v(nodes(i).number - Nx,1) - vel_v(nodes(i).number + Nx,1));
        elseif(nodes(i).code == 2)
            if(nodes(i).number == 1 )
                output_x(nodes(i).number,i) = (vel_u(nodes(i).number + 1,1));
                output_y(nodes(i).number,i) = -(vel_v(nodes(i).number + Nx,1));
            elseif(nodes(i).number == Nx)
                output_x(nodes(i).number,i) = (-vel_u(nodes(i).number - 1,1));
                output_y(nodes(i).number,i) = -(vel_v(nodes(i).number + Nx,1));
            else
                output_x(nodes(i).number,i) = (vel_u(nodes(i).number + 1,1) - vel_u(nodes(i).number - 1,1));
                output_y(nodes(i).number,i) = -(vel_v(nodes(i).number + Nx,1));
            end
        
        elseif(nodes(i).code == 5)
            if(nodes(i).number == Nx*Ny - Nx + 1)
                output_x(nodes(i).number,i) = (vel_u(nodes(i).number + 1,1));
                output_y(nodes(i).number,i) = (vel_v(nodes(i).number - Nx,1));
            elseif(nodes(i).number == Nx*Ny)
                output_x(nodes(i).number,i) = (-vel_u(nodes(i).number - 1,1));
                output_y(nodes(i).number,i) = (vel_v(nodes(i).number - Nx,1));
            else
                output_x(nodes(i).number,i) = (vel_u(nodes(i).number + 1,1) - vel_u(nodes(i).number - 1,1));
                output_y(nodes(i).number,i) = (vel_v(nodes(i).number - Nx,1));
            end
        elseif(nodes(i).code == 4)
            output_x(nodes(i).number,i) = (vel_u(nodes(i).number + 1,1));
            output_y(nodes(i).number,i) = (vel_v(nodes(i).number - Nx,1) - vel_v(nodes(i).number + Nx,1));
        elseif(nodes(i).code == 3)
            output_x(nodes(i).number,i) = (- vel_u(nodes(i).number - 1,1));
            output_y(nodes(i).number,i) = (vel_v(nodes(i).number - Nx,1) - vel_v(nodes(i).number + Nx,1));
        end
    end
    output = (1/(2*dx))*output_x + (1/(2*dy))*output_y;
    output = sparse(output);
end