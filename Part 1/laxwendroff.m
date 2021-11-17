function solnext = laxwendroff(sol,dt,dx,u0,Nx)
    a = 1;
    sigma = a*dt/dx;
    alpha = 0.5*sigma*(sigma+1);
    beta = 1-sigma^2;
    gamma = 0.5*sigma*(sigma-1);
    ui = sol(1:end-1,end);
    if(size(sol,2) >2)
        uNm1 = sol(end-1,end);
        uNm2 = sol(end-2,end);
        uN = extrapolate(uNm1, uNm2);
    else
        uN = 0;
    end
    global A
    if isempty(A)
        A1=circshift(eye(Nx-1),1)*alpha;
        A2=eye(Nx-1)*beta;
        A3=circshift(eye(Nx-1),-1)*gamma;
        A=sparse(A1+A2+A3);
        A(1,end) = 0;
        A(end,1) = 0;
    end

    
    b = [u0*alpha;zeros(Nx-3,1);uN*beta];
    unext = [A*ui+b;uN];
    solnext = [sol,unext];
end