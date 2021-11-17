clear all
close all
clc
set(0,'defaultTextInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');

Nx = 1000; % number of spatial grid points 
v = 1; % velocity of fluid
alpha = 0.5; % heat exchange parameter
L = 5; % length of the pipe
Tcool = 50; % temperature of the cooling
Thot = 200; % temperature of the hot whatever
dx = L/(Nx-1); % spatial steplength
CFL = 1; % CFL number
dt = dx*CFL; % time step length
tend = 6; % maximum time
uini = zeros(Nx,1);
uini(:) = Tcool; % specifying initial conditions
sol = uini;
%%%%%%%%% upwind
lambda = dt/dx; % CFL number
alphajminus1 = v*lambda; % factor of u_{j-1}^n
betaj = 1-lambda*v-dt*alpha; %factor of u_{j}^n
A = zeros(Nx-1); % Matrix A that is used to calculate the next timestep
A(1,1) = betaj;
for k = 2:Nx-1
    A(k,k-1) = alphajminus1;
    A(k,k) = betaj;
end
A = sparse(A);
for t = 0:dt:tend
    u0 = boundary(t);
    sol = upwind2(sol,dt,dx,u0,Nx,A);
end
%%%%%%%%%
%%%%% plot
x = linspace(0,5,length(sol(:,1)));
t = linspace(0,6,length(sol(1,:)));
[X,T] = meshgrid(t,x);
mesh(X,T,sol)
grid on
xlabel('time')
ylabel('x')
zlabel('Temperature')
%%%%% end of plot
% function that generates boundary conditions u0(t) at x = 0
function u0 = boundary(t)
Tcool = 50;
Thot = 200;
if t < 0.125
    u0 = Tcool + (Thot - Tcool) * sin(4*pi*t);
elseif 0.125 <= t && t <= 1
    u0 = Thot;
elseif t > 1
    u0 = Thot + Tcool*sin(5*pi*(t-1));
end
end