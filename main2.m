clear all
close all
clc

v = 1; % velocity of fluid
alpha = 0.5; % heat exchange parameter
L = 5; % length of the pipe
Tcool = 50; % temperature of the cooling
Thot = 200; % temperature of the hot whatever
dx = L/100; % spatial steplength
CFL = 1; % CFL number
dt = dx*CFL; % time step length
tend = 6; % maximum time
uini = zeros(100,1);
uini(:) = Tcool; % specifying initial conditions
sol = uini;
for t = 0:dt:tend
    u0 = boundary(t);
    sol = upwind2(sol,dt,dx,u0);
end
%%%%% plot
x = linspace(0,5,length(sol(:,1)));
t = linspace(0,6,length(sol(1,:)));
plot(t,sol(1,:),'bX')
grid on
xlabel('time \rightarrow')
ylabel('solution \rightarrow')
legend('solution at x = 0')
%%%%% end of plot
t = linspace(0,6,length(sol(1,:)));
for k=1:length(t)
u0test(k) = boundary(t(k));
end
% function that generates boundary conditions u0 at x = 0
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