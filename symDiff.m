clear
clc
syms f(theta,w0,w_prime,beta) F(theta, beta, A , B)
f(theta,w0,w_prime,beta) = exp(-beta*theta)*(w0*cos(beta*theta) + (w_prime/beta+w0)*sin(beta*theta));
F (theta , beta , A,B) = exp(-beta*theta)*(A*cos(beta*theta) + B*sin(beta*theta))
dr = diff(f , theta)
ddr = diff(diff(f, theta),theta);
subs(ddr , theta , 0)