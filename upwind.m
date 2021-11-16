function solnext = upwind(sol,dt,dx,u0,Nx)
lastrow = sol(:,end); %obtaining the last time step to calculate the next time step 
lambda = dt/dx; % CFL number
a=1; % transport velocity
alpha = a*dt/dx; % factor of u_{j-1}^n
beta = 1-a*dt/dx; %factor of u_{j}^n
A = zeros(Nx-1); % Matrix A that is used to calculate the next timestep
A(1,1) = beta;
for k = 2:Nx-1
    A(k,k-1) = alpha;
    A(k,k) = beta;
end
A = sparse(A);
b = zeros(Nx-1,1);
b(1) = alpha * lastrow(1); % used to get u_1 with the help of u_0 
u_next_timestep_2_Nx = A * lastrow(2:Nx)+b;
u_next_timestep = [u0;u_next_timestep_2_Nx];
solnext = [sol,u_next_timestep];
end