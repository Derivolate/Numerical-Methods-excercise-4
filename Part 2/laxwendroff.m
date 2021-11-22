function solnext = laxwendroff(sol,u0)
    global dt; 
    global dx;
    global Nx;
    global alpha
    global v;
    global Tcool;
    global A;
    
    a = v;
    b = alpha;
    c = Tcool;
    sigma = a*dt/dx;
    
    chi = sigma/2*(sigma+(1-b*dt));
    phi = sigma/2*(sigma-(1-b*dt));
    theta = c*dt*(1-b*dt/2)/2;
    ui = sol(1:end-1,end);
    if(size(sol,2) >2)
        uNm1 = sol(end-1,end);
        uNm2 = sol(end-2,end);
        uN =2*uNm1 - uNm2;
    else
        uN = 0;
    end
    if isempty(A)
        psi = 1-b*dt*(1-b*dt/2)-sigma^2;
        A1=circshift(eye(Nx-1),1)*chi;
        A2=eye(Nx-1)*psi;
        A3=circshift(eye(Nx-1),-1)*phi;
        A=sparse(A1+A2+A3);
        A(1,end) = 0;
        A(end,1) = 0;
    end
    b = [u0*chi+theta;ones(Nx-3,1)*theta;uN*phi+theta];
    A = sparse(A);
    size(ui)
    size(A)
    size(b)
    unext = [A*ui+b;uN];
    solnext = [sol,unext];
end
