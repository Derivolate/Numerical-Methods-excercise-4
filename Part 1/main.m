clear all
close all
clc
set(0,'defaultTextInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
D = 4.5;
a = 1;
tau = 2.2;
Nx = 100; % number of spatial grid points
dx = D/(Nx-1);
CFL = 0.5;
dt = dx*CFL;
tend = 6;
uini = zeros(Nx,1);
global A
%%% calling upwind sin
clear global A
global A
sol_upwind_sin = uini;
for t = dt:dt:tend
    u0 = gsin(t,tau);
    sol_upwind_sin = upwind(sol_upwind_sin,dt,dx,u0,Nx);
end
%%% calling Lax Friedrich sin
clear global A
global A
sol_laxfriedrich_sin = uini;
for t = dt:dt:tend
    u0 = gsin(t,tau);
    sol_laxfriedrich_sin = laxfriedrich(sol_laxfriedrich_sin,dt,dx,u0,Nx);
end
%%% calling Lax Wendroff sin
clear global A
global A
sol_laxwendroff_sin = uini;
for t = dt:dt:tend
    u0 = gsin(t,tau);
    sol_laxwendroff_sin = laxwendroff(sol_laxwendroff_sin,dt,dx,u0,Nx);
end
%%%% implement exact solution sin
x = linspace(0,D,Nx);
t = linspace(0,tend,tend/(CFL*dx)+1);
uexactsin = zeros(length(x),length(t));
for jj = 1:length(t)
    for kk = 1:Nx
        uexactsin(kk,jj) = gsin(t(jj)-a*x(kk),tau);
    end
end
%%%%
%%% calling upwind sq
clear global A
global A
sol_upwind_sq = uini;
for t = dt:dt:tend
    u0 = gsq(t,tau);
    sol_upwind_sq = upwind(sol_upwind_sq,dt,dx,u0,Nx);
end
%%% calling Lax Friedrich sq
clear global A
global A
sol_laxfriedrich_sq = uini;
for t = dt:dt:tend
    u0 = gsq(t,tau);
    sol_laxfriedrich_sq = laxfriedrich(sol_laxfriedrich_sq,dt,dx,u0,Nx);
end
%%% calling Lax 
% Wendroff sq
clear global A
global A
sol_laxwendroff_sq = uini;
for t = dt:dt:tend
    u0 = gsq(t,tau);
    sol_laxwendroff_sq = laxwendroff(sol_laxwendroff_sq,dt,dx,u0,Nx);
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

%%%%% plot
x = linspace(0,D,Nx);
hold on
subplot(2,1,1)
plot(x,uexactsin(:,end),'k',x,sol_upwind_sin(:,end),'b',x,sol_laxfriedrich_sin(:,end),'c',x,sol_laxwendroff_sin(:,end),'r')
legend('exact solution', 'upwind','Lax-Friedrich','Lax-Wendroff')
grid on
xlabel('x')
ylabel('u')
subplot(2,1,2)
plot(x,uexactsq(:,end),'k',x,sol_upwind_sq(:,end),'b',x,sol_laxfriedrich_sq(:,end),'c',x,sol_laxwendroff_sq(:,end),'r')
legend('exact solution', 'upwind','Lax-Friedrich','Lax-Wendroff')
grid on
xlabel('x')
ylabel('u')
%%%%% end of plot
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