dy = -1;
theta = (0:0.001:2*pi)';
x = cos(theta);
y = sin(theta);
dr_dtheta = (x + y*dy)./(dy*x - y);
close all
daspect([1 , 1,  1])
plot(rad2deg(theta) , dr_dtheta)