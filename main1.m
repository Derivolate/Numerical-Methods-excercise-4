clear all
close all
clc

D = 4.5;
a = 1;
tau = 2.2;
dx = D/100;
CFL = 1;
dt = dx*CFL;
tend
uini = zeros(100,1);
sol = uini;

for t = 0:dt:tend
    u0 = gsin(t,tau);
    sol = upwind(sol,dt,dx,u0);
end

%Initial condition functions
function u0 = gsin(t, tau)
    u0 = sin((t/tau)*2*pi);
end
function u0 = gsq(t,tau) 
    u0 = square((t/tau)*2*pi);
end