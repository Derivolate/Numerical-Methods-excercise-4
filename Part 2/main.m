clear all
close all
clc
set(0,'defaultTextInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');

global Nx;
global v;
global alpha; 
global L;
global Tcool;
global Thot;
global dx;
global dt;

Nx = 100; % number of spatial grid points 
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
sol_upwind = uini;
sol_laxwendroff = uini;
%%%%%%%%% upwind
global A
clear A;
for t = dt:dt:tend
    u0 = boundary(t);
    sol_upwind = upwind(sol_upwind,u0);
    %sol_laxwendroff = laxwendroff(sol_laxwendroff,u0);
end

%%%%%%%%%
%%%%% plot for (a)
% x = linspace(0,5,length(sol_upwind(:,1)));
% t = linspace(0,6,length(sol_upwind(1,:)));
% [X,T] = meshgrid(t,x);
% subplot(1,2,1)
% mesh(X,T,sol_upwind)
% grid on
% xlabel('time')
% ylabel('x')
% zlabel('Temperature')
% mesh(X,T,sol_laxwendroff)
% grid on
% xlabel('time')
% ylabel('x')
% zlabel('Temperature')
% subplot(1,2,2)
%%%%% end of plot for (a)

%%%% plot for (b)
x = linspace(0,5,length(sol_upwind(:,1)));
subplot(1,2,1)
t_3 = round(length(sol_upwind(1,:))/2)
t_6 = length(sol_upwind(1,:))
plot(x,sol_upwind(:,t_3));
grid on
hold on
title('t = 3')
xlabel('position x')
ylabel('Temperature T')
legend('upwind','Lax Wendroff')
subplot(1,2,2)
plot(x,sol_upwind(:,t_6));
grid on
hold on
title('t = 6')
xlabel('position x')
ylabel('Temperature T')
legend('upwind','Lax Wendroff')

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