clear all
close all
clc

%================== USER INPUT DATA ===================%
nt = 200;            % max time steps 
Nx = 40; Ny = 40;    % Number of grids
dt = 1e-2;           % time step;
visRate = 1;         % zagęszczenie wektorów prędkości
N = 400;             % zagęszczenie linii prądu
%======================================================%


% domain size
Lx = 1; Ly = 1; 
%viscosity
viscosity = 2e-2; 
% Grid size (Equispaced)
dx = Lx/Nx;
dy = Ly/Ny;
x = (1:Nx)*dx;
y = (1:Ny)*dy;

u = zeros(Nx, Ny);
v = zeros(Nx, Ny);
U = zeros(Nx-2, Ny-2);
V = zeros(Nx-2, Ny-2);

nodes = grid_fill(Nx,Ny,Lx,Ly);
flow = struct('vel_u',{},...
    'vel_v',{},... 
    'vel_u_tylda',{},... 
    'vel_v_tylda',{},... 
    'vel_u_star',{},...
    'vel_v_star',{},...
    'pressure',{},...  
    'pressure_star',{},... 
    'pressure_tylda',{},...
    'tmp_u',{},...
    'tmp_v',{},...
    'correction_factor',{});

data = struct('res_1',{},...
    'res_2',{},...
    'res_3',{});

% Coordinate of each grid (cell corners)
xce = ((1:Nx-2))*dx;
yce = ((1:Ny-2))*dy;

cm = hsv(ceil(100/0.7));  cm = flipud(cm(1:100,:));
[Xce,Yce] = meshgrid(xce,yce); % mesh
[U_post,V_post] = postprocessor_data_conv(U,V,Nx-2,Ny-2);
figure(1);
[~,h_abs] = contourf(Xce',Yce',sqrt(U_post.^2+V_post.^2)','LineColor','none');
title('Velocity magnitude');  xlabel('x-location');  ylabel('y-location')  
axis('equal',[dx Lx-2*dx dy Ly-2*dy]);  colormap(cm);  colorbar('westoutside'); 
 
% Discretization step
laplace_domena = laplace_operator_domena(Nx,Ny,dx,dy,nodes);
[pressure_gradient_x, pressure_gradient_y] = grad_pressure(Nx,Ny,dx,dy);
laplace_ghost_nodes = Neumann_condition(Nx,Ny,Lx,Ly,nodes);
[grad_x, grad_y] = gradient_operator(Nx,Ny,dx,dy);


time = dt;
for iter = 1:nt
    [flow(iter).vel_u, flow(iter).vel_v] = boundary_conditions(Lx,Ly,Nx,Ny,u,v,time);
    if(iter < 3)
    % Initialization
    flow(iter).pressure = zeros(Nx*Ny,1);
    flow(iter).correction_factor = zeros(Nx*Ny,1);
    flow(iter).vel_u_tylda = zeros(Nx*Ny,1);
    flow(iter).vel_v_tylda = zeros(Nx*Ny,1);
    flow(iter).tmp_u = zeros(Nx*Ny,1);
    flow(iter).tmp_v = zeros(Nx*Ny,1);
    flow(iter).vel_u_star = zeros(Nx,Ny);
    flow(iter).vel_v_star = zeros(Nx,Ny);

    else
    %================= STEP 1 =================%
    flow(iter).pressure_star = flow(iter - 1).pressure;
    flow(iter).vel_u_star = flow(iter - 1).tmp_u;
    flow(iter).vel_v_star = flow(iter - 1).tmp_v;
    flow(iter).pressure_tylda = flow(iter).pressure_star + (4/3.*flow(iter - 1).correction_factor) - (1/3.*flow(iter - 2).correction_factor);
    %==========================================%
    %========== STEP 2 - x direction ==========%
    [uUx, vUy] = NS_term1(Nx,Ny,dx,dy,nodes,flow(iter).vel_u_star,flow(iter).vel_v_star);
    ns_term2 = 1/2 * NS_term2(Nx,Ny,dx,dy,flow(iter).vel_u_star,flow(iter).vel_v_star,nodes);

    RHS_step2_u = -pressure_gradient_x * flow(iter).pressure_tylda + 2/dt * flow(iter - 1).tmp_u - 1/(2*dt) * flow(iter - 2).tmp_u;
    LHS_step2_u = 3/(2*dt)*speye(Nx*Ny) - viscosity*laplace_domena + ns_term2 + uUx + vUy ;
    
    % preconditioner
    setup = struct('type','ilutp','droptol',1e-6);
    [L,U] = ilu(LHS_step2_u,setup);

    % correction of equation factors
    [LHS_step2_u,RHS_step2_u] = equation_correction_x(Lx,Ly,Nx,Ny,u,v,nodes,LHS_step2_u,RHS_step2_u,time);
    [flow(iter).vel_u_tylda, fl1, data(iter).res_1, it1, rv1] = cgs(LHS_step2_u,RHS_step2_u,1e-6,70,L,U);
    flow(iter).tmp_u = flow(iter).vel_u_tylda;
    %==========================================%
    %========== STEP 2 - y direction ==========%
    RHS_step2_v = -pressure_gradient_y * flow(iter).pressure_tylda + 2/dt * flow(iter - 1).tmp_v - 1/(2*dt) * flow(iter - 2).tmp_v;
    LHS_step2_v = 3/(2*dt)*speye(Nx*Ny) - viscosity*laplace_domena + uUx + vUy + ns_term2;
    
    % preconditioner
    [L,U] = ilu(LHS_step2_v,setup);
    
    % correction of equation factors
    [LHS_step2_v,RHS_step2_v] = equation_correction_y(Lx,Ly,Nx,Ny,u,v,nodes,LHS_step2_v,RHS_step2_v,time);
    [flow(iter).vel_v_tylda, fl2, data(iter).res_2, it2, rv2] = pcg(LHS_step2_v,RHS_step2_v,1e-6,100,L,U);
    flow(iter).tmp_v = flow(iter).vel_v_tylda;
    %==========================================%
    %============= data converting ============%
    flow(iter).vel_u_tylda = reshape(flow(iter).vel_u_tylda',Nx,Ny);
    flow(iter).vel_v_tylda = reshape(flow(iter).vel_v_tylda',Nx,Ny);
    %==========================================%
    %================= STEP 3 =================%
    RHS_step3 = 3/(2*dt) * (grad_x * flow(iter).tmp_u + grad_y * flow(iter).tmp_v); 
    LHS_step3 = laplace_ghost_nodes;

    setup = struct('type','ilutp','droptol',1e-6);
    [L,U] = ilu(LHS_step3, setup);
    [flow(iter).correction_factor, fl3, data(iter).res_3, it3, rv3] = cgs(LHS_step3,RHS_step3,1e-6,100, L, U);
    %==========================================%
    %================= STEP 4 =================%
    flow(iter).pressure = flow(iter).pressure_star + flow(iter).correction_factor - viscosity * (grad_x * flow(iter).tmp_u + grad_y * flow(iter).tmp_v);
    pressure = reshape(flow(iter).pressure',Nx,Ny);
    %==========================================%
    %============= data converting ============%
    for i = 1:Nx-2
        U_post(i,:) = flow(iter).vel_u_tylda(i+1,2:end-1);
        V_post(i,:) = flow(iter).vel_v_tylda(i+1,2:end-1);
    end
    %==========================================%
    
    %============= CONTOUR UPDATING ===========%
        [U_post,V_post] = postprocessor_data_conv(U_post,V_post,Nx-2,Ny-2);
        h_abs.ZData = sqrt(V_post.^2 + U_post.^2);
        drawnow
        title('time = ',num2str(time));
    %==========================================%
    end
    time = time + dt;
end


%=========================== POSTPROCESSING ===========================%
% Downsample the data
figure(6);
xced = xce(1:visRate:end);
yced = yce(1:visRate:end);
[Xced,Yced] = meshgrid(xced, yced);

uced = U_post(1:visRate:end,1:visRate:end);
vced = V_post(1:visRate:end,1:visRate:end);
h_quiver = quiver(Xced',Yced',uced,-vced,3,'Color',[0,0,0]);
harrow = annotation('textarrow',[0.3 0.7],[0.96 0.96],"LineWidth",2);


cm = hsv(ceil(100/0.7));  cm = flipud(cm(1:100,:));
figure(2);  contourf(xce,yce,U_post',23,'LineColor','none');
title('U-velocity');  xlabel('x-location');  ylabel('y-location')  
axis('equal',[dx Lx-2*dx dy Ly-2*dy]);  colormap(cm);  colorbar('westoutside'); 

figure(3);  contourf(xce,yce,V_post',23,'LineColor','none');
title('V-velocity');  xlabel('x-location');  ylabel('y-location')  
axis('equal',[dx Lx-2*dx dy Ly-2*dy]);  colormap(cm);  colorbar('westoutside'); 

cm = hsv(ceil(100/0.7));  cm = flipud(cm(1:100,:));
figure(4);  contourf(xce,yce,sqrt(U_post.^2+V_post.^2)',23,'LineColor','none');
title('Velocity magnitude');  xlabel('x-location');  ylabel('y-location')  
axis('equal',[dx Lx-2*dx dy Ly-2*dy]);  colormap(cm);  colorbar('westoutside'); 

xstart = max(xce)*rand(N,1);  ystart = max(y)*rand(N,1);
[X,Y] = meshgrid(xce,yce);
figure(5);  h=streamline(X,Y,U_post',-V_post',xstart,ystart,[0.1, 200]);
title('Streamlines');  xlabel('x-location');  ylabel('y-location')
axis('equal',[0 Lx 0 Ly]);  set(h,'color','k')

pressure = postprocessor_data_conv(pressure,pressure,Nx,Ny);
figure(7);  contourf(x,y,pressure',23,'LineColor','none');
title('Pressure');  xlabel('x-location');  ylabel('y-location')
axis('equal',[dx Lx-2*dx dy Ly-2*dy]);  colormap(cm);  colorbar('westoutside'); 

res_1 = zeros(iter,1);
res_2 = zeros(iter,1);
res_3 = zeros(iter,1);

for i = 3:iter
    res_1(i,1) = data(i).res_1;
    res_2(i,1) = data(i).res_2;
    res_3(i,1) = data(i).res_3;
end

figure(8)
semilogy(1:iter,res_1,'-o')
hold on
semilogy(1:iter,res_2,'-o')
hold on
semilogy(1:iter,res_3,'-o')
legend('Diffusion equation - direction x','Diffusion equation - direction y','Projection step','Location','southeast')
xlabel('Iteration number')
ylabel('Relative residual')






