clear
clc

tyre_sim_structs
%%



%% plots
close all
plot(undeformed_tyre.x , undeformed_tyre.y , "--")
hold on 
for beta = linspace(0.1 , 10 , 10)
deformed_tyre = undeformed_tyre;
deformed_tyre = getDeformedRadius(penetration_end_ind, 1,0.1, 0, deformed_tyre);
deformed_tyre = getDeformedRadius(penetration_start_ind , -1 , 0.1 , 0 , deformed_tyre);
deformed_tyre.x = deformed_tyre.radius .* cos(theta);
deformed_tyre.y = deformed_tyre.radius .* sin(theta);
plotDeformedTyre(deformed_tyre, gcf(), penetration_start_ind, penetration_end_ind)
end
daspect([1 1 1])
legend 
%% functions
% assuming end is always after start, meaning the contact patch goes CCW
% from start to end
function tyre = getDeformedRadius(...
                                                                                penetration_ind, ...
                                                                                direction, ...
                                                                                penetration, ...
                                                                                penteration_angle , ...
                                                                                tyre)
global beta
num_nodes = length(tyre.radius);
deformed_radius = zeros(1 ,num_nodes);
current_ind = penetration_ind;

 = 2*pi/num_nodes;
current_theta = 0;
while abs(current_ind - penetration_ind) < 100
    tyre.radius(mod(current_ind , num_nodes) + 1) = tyre.radius(mod(current_ind , num_nodes) +1) - ...
        exp(-beta * current_theta) * (penetration* cos(beta*current_theta) + (penteration_angle/beta  - penetration) * sin(beta * current_theta));
    fprintf("deformed_radius (%i) = %.3f\n", mod(current_ind , num_nodes) + 1, deformed_radius(mod(current_ind , num_nodes)+1));
    current_theta = current_theta + dTheta;
    current_ind = current_ind + direction;
end    
end
function plotDeformedTyre(tyre, fig_handle, penetration_start_ind, penetration_end_ind)
    figure(fig_handle);
    hold on 
    plot(tyre.x(penetration_end_ind+1:end), tyre.y(penetration_end_ind+1:end));
    plot(tyre.x(1:penetration_start_ind), tyre.y(1:penetration_start_ind));
    plot(tyre.x([penetration_start_ind , penetration_end_ind+1]), tyre.y([penetration_start_ind, penetration_end_ind+1]));
end

