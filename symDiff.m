clear
clc
syms f(theta,w0,w_prime,beta) F(theta, beta, A , B)
f(theta,w0,w_prime,beta) = exp(-beta*theta)*(w0*cos(beta*theta) + (w_prime/beta+w0)*sin(beta*theta));
F (theta , beta , A,B) = exp(-beta*theta)*(A*cos(beta*theta) + B*sin(beta*theta))
dr = diff(f , theta)
ddr = diff(diff(f, theta),theta);
subs(ddr , theta , 0)
%%
clear
clc
close all
load poly2sol.mat
syms x1 x2 x3 y1 y2 y3 dy3
X1 = 0;
Y1 = 0;
X2 = 5;
Y2 = -5;
X3 = 1;
Y3 = 5;
get_val = @(a)eval(subs(a , [x1,x2,x3, y1,y2,y3,dy3], [X1, X2, X3, Y1,Y2,Y3, 0]));
p1_x = linspace(X1, X3, 100);
p2_x = linspace(X3 , X2, 100);
p1_y = polyval([get_val(sol.a1), get_val(sol.b1), get_val(sol.c1)], p1_x);
p2_y = polyval([get_val(sol.a2), get_val(sol.b2), get_val(sol.c2)], p2_x);
plot(p1_x, p1_y);
hold on 
plot(p2_x , p2_y);
%%
clear all
clc
syms a  b xl xr yl yr y0 x
A = [xl^2, xl;xr^2, xr]
Y = [yl-y0;yr-y0]
C = A\Y
%%
clear
clc
syms x x0 beta A B
beam = exp(-beta*x)*(A * cos(beta*x) + B*sin(beta*x));
beam_offset = expand(subs(beam , x + x0));
pretty(collect(beam_offset, exp(-beta*x)*exp(-beta*x0)))
