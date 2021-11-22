%% plotcalling for b
clc
clear all
dx = input('dx: ');
CFL = input('CFL: ');
fprintf(['Different methods' newline '1 Upwind' newline '2 Lax-Friedrich' newline '3 Lax-Wendroff' newline])
method = input('Method:');
fprintf(['Different boundary conditions' newline '1 sin' newline '2 square'])
bound = input('Boundary condition:');

AA = mainb(dx,CFL,method,bound);
function doplot = mainb(dx,CFL,method,bound)

D = 4.5;
a = 1;
tau = 2.2;
Nx = 1+D/dx; % number of spatial grid points
dt = dx*CFL;
tend = 6;
uini = zeros(Nx,1);
grid on
xlabel('x')
ylabel('u')
global A
clear global A
    global A
if method == 1 && bound == 1
    %%% calling upwind sin
    
    sol_upwind_sin = uini;
    for t = dt:dt:tend
        u0 = gsin(t,tau);
        sol_upwind_sin = upwind(sol_upwind_sin,dt,dx,u0,Nx);
    end
    %%%%% plot
    x = linspace(0,D,Nx);
    hold on
    plot(x,sol_upwind_sin(:,end))
elseif method == 2 && bound == 1
    %%% calling Lax Friedrich sin
    sol_laxfriedrich_sin = uini;
    for t = dt:dt:tend
        u0 = gsin(t,tau);
        sol_laxfriedrich_sin = laxfriedrich(sol_laxfriedrich_sin,dt,dx,u0,Nx);
    end
    %%%%% plot
    x = linspace(0,D,Nx);
    hold on
    plot(x,sol_laxfriedrich_sin(:,end))
elseif method == 3 && bound == 1
    %%% calling Lax Wendroff sin
    sol_laxwendroff_sin = uini;
    for t = dt:dt:tend
        u0 = gsin(t,tau);
        sol_laxwendroff_sin = laxwendroff(sol_laxwendroff_sin,dt,dx,u0,Nx);
    end
    %%%%% plot
    x = linspace(0,D,Nx);
    hold on
    plot(x,sol_laxwendroff_sin(:,end))
elseif method == 1 && bound == 2
    %%% calling upwind sq
    clear global A
    global A
    sol_upwind_sq = uini;
    for t = dt:dt:tend
        u0 = gsq(t,tau);
        sol_upwind_sq = upwind(sol_upwind_sq,dt,dx,u0,Nx);
    end
    %%%%% plot
    x = linspace(0,D,Nx);
    hold on
    plot(x,sol_upwind_sq(:,end))
elseif method == 2 &&  bound == 2
    %%% calling Lax Friedrich sq
    sol_laxfriedrich_sq = uini;
    for t = dt:dt:tend
        u0 = gsq(t,tau);
        sol_laxfriedrich_sq = laxfriedrich(sol_laxfriedrich_sq,dt,dx,u0,Nx);
    end
    x = linspace(0,D,Nx);
    hold on
    plot(x,sol_laxfriedrich_sq(:,end))

elseif method == 3 &&  bound == 2
    %%% calling Lax
    % Wendroff sq
    sol_laxwendroff_sq = uini;
    for t = dt:dt:tend
        u0 = gsq(t,tau);
        sol_laxwendroff_sq = laxwendroff(sol_laxwendroff_sq,dt,dx,u0,Nx);
    end
    x = linspace(0,D,Nx);
    hold on
    plot(x,sol_laxwendroff_sq(:,end))
else
    error('Wrong input!')
end

%%%% implement exact solution sq
x = linspace(0,D,Nx);
t = linspace(0,tend,tend/(CFL*dx)+1);
uexactsq = zeros(length(x),length(t));
for jj = 1:length(t)
    for kk = 1:Nx
        uexactsq(kk,jj) = gsq(t(jj)-a*x(kk),tau);
    end
end
%%%%
%Initial condition functions
    function u0 = gsin(t, tau)
        if t > 0
            u0 = sin((t/tau)*2*pi);
        elseif t <= 0
            u0 = 0;
        end
    end

    function u0 = gsq(t,tau)
        if t > 0
            u0 = square((t/tau)*2*pi);
        elseif t <= 0
            u0 = 0;
        end
    end
doplot = 1;
end