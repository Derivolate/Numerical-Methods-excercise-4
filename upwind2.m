function solnext = upwind2(sol,dt,dx,u0,Nx)
lastrow = sol(:,end); %obtaining the last time step to calculate the next time step
lambda = dt/dx; % CFL number
v=1; % velocity of the fluid
alpha = 0.5; % heat exchange parameter, note that alpha =/= alphaj-1 !
Tcool = 50;
alphajminus1 = v*lambda; % factor of u_{j-1}^n
betaj = 1-lambda*v-dt*alpha; %factor of u_{j}^n
A = zeros(Nx-1); % Matrix A that is used to calculate the next timestep
A(1,1) = betaj;
for k = 2:Nx-1
    A(k,k-1) = alphajminus1;
    A(k,k) = betaj;
end
A = sparse(A);
b = zeros(Nx-1,1);
b(2:Nx-1,1) = dt*alpha*Tcool;
b(1,1) = alphajminus1*lastrow(1) + dt*alpha*Tcool;
u_next_timestep_2_Nx = A * lastrow(2:Nx)+b;
u_next_timestep = [u0;u_next_timestep_2_Nx];
solnext = [sol,u_next_timestep];
end