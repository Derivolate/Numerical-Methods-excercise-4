function solnext = upwind(sol,u0,A)
    global dt; 
    global dx;
    global Nx;
    global alpha
    global v;
    global Tcool;
    global A;
    lastrow = sol(:,end); %obtaining the last time step to calculate the next time step
    lambda = dt/dx; % CFL number
    alphajminus1 = v*lambda; % factor of u_{j-1}^n
    betaj = 1-lambda*v-dt*alpha; %factor of u_{j}^n
    b = zeros(Nx-1,1);
    b(2:Nx-1,1) = dt*alpha*Tcool;
    b(1,1) = alphajminus1*lastrow(1) + dt*alpha*Tcool;
    
    %If A does not exist in the global namespace, define it
    if isempty(A)
        A = sparse(zeros(Nx-1)); % Matrix A that is used to calculate the next timestep
        A(1,1) = betaj;
        for k = 2:Nx-1
            A(k,k-1) = alphajminus1;
            A(k,k) = betaj;
        end
    end
    u_next_timestep_2_Nx = A * lastrow(2:Nx)+b;
    u_next_timestep = [u0;u_next_timestep_2_Nx];
    solnext = [sol,u_next_timestep];
    
    
end