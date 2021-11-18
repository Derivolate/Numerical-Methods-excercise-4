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