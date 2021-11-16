function solnext = up(sol,dt,dx,u0)
lastrow = sol(:,end); %obtaining the last time step to calculate the next time step 
lambda = dt/dx; % CFL number
a=1; % transport velocity
alpha = a*dt/dx; % factor of u_{j-1}^n
beta = 1-a*dt/dx; %factor of u_{j}^n
A = zeros(99); % Matrix A that is used to calculate the next timestep
A(1,1) = beta;
for k = 2:99
    A(k,k-1) = alpha;
    A(k,k) = beta;
end
A = sparse(A);
b=zeros(99,1);
b(1) = alpha * lastrow(1);
u_next_timestep_2_100 = A * lastrow(2:100)+b;
u_next_timestep = [u0;u_next_timestep_2_100];
solnext = [sol,u_next_timestep];
end