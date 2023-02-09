function [output_u, output_v] = postprocessor_data_conv(U,V,Nx,Ny)
    for iter = 1:Ny
        output_u(:,iter) = U(:,(Ny+1)-iter); 
        output_v(:,iter) = V(:,(Ny+1)-iter); 
    end
end

