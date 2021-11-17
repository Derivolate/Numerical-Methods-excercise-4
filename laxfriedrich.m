function solnext = laxfriedrich(sol,dt,dx,u0,Nx)
    a = 1;
    sigma = a*dt/dx;
    alpha = 0.5*(1+sigma);
    beta = 0.5*(1-sigma);
    ui = sol(2:end-1,end);
    if(size(sol,2) >2)
        uNm1 = sol(end,end-1);
        uNm2 = sol(end,end-2);
        uN = extrapolate(uNm1, uNm2);
    else
        uN = 0;
    end
    A1=circshift(eye(Nx-2),1)*alpha;
    A2=circshift(eye(Nx-2),-1)*beta;
    A=(A1+A2);
    A(1,end) = 0;
    A(end,1) = 0;
    b = [u0*alpha;zeros(Nx-4,1);uN*beta];
    A = sparse(A);
    unext = [u0;A*ui+b;uN];
    solnext = [sol,unext];
end