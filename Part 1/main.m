clear all
close all
clc
D = 4.5;
a = 1;
tau = 2.2;

Nx = 1000; % number of spatial grid points
dx = D/(Nx-1);
CFL = 0.95;
dt = dx*CFL;
tend = 6;
%uini = zeros(Nx,1);
uini = zeros(Nx,1);
sol = uini;
%%% calling the upwind method
global A
clear A;
for t = dt:dt:tend
    u0 = gsq(t,tau);
    sol = laxwendroff(sol,dt,dx,u0,Nx);
end
%%%% implement exact solution
x = linspace(0,D,Nx);
t = linspace(0,tend,tend/D*Nx);
uexact = zeros(length(x),length(t));
for jj = 1:length(t)
    for kk = 1:Nx
        uexact(kk,jj) = gsin(t(jj)-a*x(kk),tau);
    end
end
%%%%
%%%%% plot
x = linspace(0,D,Nx);
plot(x,uexact(:,end),'k')
hold on
plot(x,sol(:,end))
grid on
xlabel('x \rightarrow')
ylabel('solution \rightarrow')
legend('exact solution', 'upwind')
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