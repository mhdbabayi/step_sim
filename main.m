clear
clc
close all
tyre_sim_structs
background_plots = [];
frame_plots = [];
road_s = createRoad(road_s, constants_s);
main_fig= figure();
subplot(2 , 1 , 1)
xlim([road_s.x(1) , road_s.x(end)]);
ylim([tyre_s.y_centre - tyre_s.free_radius*1.1 , tyre_s.y_centre + tyre_s.free_radius*1.2]);
hold on
daspect([1 1 1])
background_plots = drawRoad(road_s, background_plots);
% vidRecorder = VideoWriter("beamVid");
% vidRecorder.open();
while tyre_s.x_centre < (constants_s.step_position + constants_s.tyre_radius)
    try
        frame_plots = drawUndeformedTyre(tyre_s, frame_plots);
        [tyre_s, road_s] = tyreRoadContact(tyre_s , road_s , constants_s);
        frame_plots = drawPenetrations(tyre_s, road_s , frame_plots);
        direction = -1;
        r = getDeformedTyre(tyre_s , road_s, 1 ,direction , constants_s);
        plot_range = mod([tyre_s.contact_boundary_inds(1 , 1):direction:tyre_s.contact_boundary_inds(1 , 1) + 100*direction]-1, 360) + 1;
        frame_plots(end+1) = plot(tyre_s.x_centre + r(plot_range).*cos(tyre_s.theta(plot_range)) , tyre_s.y_centre+ r(plot_range).*sin(tyre_s.theta(plot_range)),LineWidth=3);
        direction = 1;
        r = getDeformedTyre(tyre_s , road_s, 1 ,direction , constants_s);
        plot_range = mod([tyre_s.contact_boundary_inds(1 , 2):direction:tyre_s.contact_boundary_inds(1 , 2) + 100*direction] -1, 360) + 1;
        frame_plots(end+1) = plot(tyre_s.x_centre + r(plot_range).*cos(tyre_s.theta(plot_range)) , tyre_s.y_centre+ r(plot_range).*sin(tyre_s.theta(plot_range)),LineWidth=3);
%         vidRecorder.writeVideo(getframe);
        pause();
        delete(frame_plots);
        frame_plots = [];
    catch
    end
    tyre_s.x_centre = tyre_s.x_centre + 0.01;
end
% vidRecorder.close();
clear road tyre constants
%% functions
% main function
function [tyre , road] = tyreRoadContact(tyre ,road, constants)
[tyre, road] = getAllPenetrations(tyre , road);
[tyre, road] = getContactCentres(tyre, road);
[tyre ,road] = getContactBoundaries(tyre, road, constants);
end

% helper functions, functions that don't return objects
function is_in_circle = getCirclePenetration(centre_x , centre_y, radius, point_x , point_y)
is_in_circle = ((point_x - centre_x)^2 + (point_y - centre_y)^2) < radius^2;
end

function first_derivative = getCentralDiffDerivative(x , y)
first_derivative = (y(3:end) - y(1:end-2))./(x(3:end) - x(1:end-2));
first_derivative = [first_derivative(1);first_derivative;first_derivative(end)];
end

function penetration_inds = findPenetratingInds(x_centre, y_centre,x_road,y_road, radius, search_start, search_end)
% searches in the range search_start to search_end and gets the boundaries
% of the first penetration
search_ind = search_start;
penetration_inds = [];
while search_ind < search_end
    if getCirclePenetration(x_centre , y_centre, radius , x_road(search_ind) , y_road(search_ind))
        penetration_inds(1,1)= search_ind;
        while getCirclePenetration(x_centre , y_centre, radius , x_road(search_ind) , y_road(search_ind))
            search_ind = search_ind + 1;
        end
        penetration_inds(1 , 2) = search_ind;
        return
    end
    search_ind = search_ind + 1;
end
end

function tyre_ind = getTyreIndexfromRoadIndex(tyre , road, road_ind)
x = road.x(road_ind) - tyre.x_centre;
y = road.y(road_ind) - tyre.y_centre;
theta = mod(atan2(y , x) , 2*pi);
tyre_ind = interp1(tyre.theta , [1:length(tyre.theta)], theta, "nearest");
end

function polar_coordinates = road2tyreTF(tyre, road_point_x, road_point_y)
polar_coordinates(1 ,1) = sqrt((road_point_x - tyre.x_centre)^2 + (road_point_y - tyre.y_centre)^2);
polar_coordinates(2 ,1) = mod(atan2(road_point_y - tyre.y_centre , road_point_x - tyre.x_centre), 2*pi);
end

function polar_derivative = road2tyreDerivativeTF(tyre, road, road_point_ind, direction)
x = road.x(road_point_ind) -tyre.x_centre;
y = road.y(road_point_ind) - tyre.y_centre;
R = sqrt(x^2 + y^2);
theta = atan2(y , x);
y_prime = road.gradient(road_point_ind);
polar_derivative = direction * (x + y*y_prime)/(y_prime * cos(theta) - sin(theta));
%     polar_derivative = (1/R)*(x^2 - y^2 * y_prime^2 )/(y_prime);
%     polar_derivative = (x + y*y_prime)/(x*y_prime - y);
end


function polar_derivative = polarDerivative(x , y , dy)
R = sqrt(x^2 + y^2);
polar_derivative = R*(x + y*dy)/(x*dy - y);
end

function polar_second_derivative  = polarSecondDerivative(x , y , dy , ddy)
r = sqrt(x^2 + y^2);
theta = atan2(y , x);
dr_dtheta = polarDerivative(x , y , dy);
dx_dtheta = dr_dtheta*x/r  - y;
dy_dtheta = dr_dtheta*y/r + x;
dydx_dtheta = ddy * dx_dtheta + (ddy/dy)*dy_dtheta; % d(dy/dx)/dtheta;
polar_second_derivative = ((x*dy - y)*(dx_dtheta + y*dydx_dtheta + dy*dy_dtheta) - ...
    (y*dy + x)*(dy*dx_dtheta + dydx_dtheta*x - dy_dtheta)) * r/ ...
    (x*dy - y)^2 +  dr_dtheta*(y*dy  +x)/(x*dy - y);

end

function global_coordinates = tyre2roadTF(tyre, theta , r)
global_coordinates(1 , 1) = tyre.x_centre + r*cos(theta);
global_coordinates(2 , 1) = tyre.y_centre + r*sin(theta);
end

function deformed_tyre = getDeformedTyre(tyre, road, contact_ind, direction, constants)
deformed_tyre = tyre.r_free;
% direction 1 means fore and direction -1 means aft
if direction == 1
    current_ind_tyre = tyre.contact_boundary_inds(contact_ind , 2);
    current_ind_road = road.contact_boundary_inds(contact_ind , 2);
elseif direction == -1
    current_ind_tyre = tyre.contact_boundary_inds(contact_ind , 1);
    current_ind_road = road.contact_boundary_inds(contact_ind , 1);
else
    error("direction should 1 or -1");
end
num_nodes = length(tyre.theta);
dTheta = 2*pi/num_nodes;
current_theta = dTheta;
initial_point_polar = road2tyreTF(tyre, road.x(current_ind_road), ...
    road.y(current_ind_road));
initial_penetration =  tyre.free_radius - initial_point_polar(1) ;
%     initial_theta =  pi/2 - initial_point_polar(2) + direction * (atan(road.gradient(current_ind_road)));
initial_theta = -road2tyreDerivativeTF(tyre, road, current_ind_road, direction);
deformed_tyre(current_ind_tyre) = initial_point_polar(1);
current_ind_tyre = current_ind_tyre + direction;
while(exp(-constants.beta*current_theta) > 0.01)
    deformed_tyre(mod(current_ind_tyre -1 , num_nodes) + 1) = tyre.free_radius -...
        exp(-constants.beta * current_theta)*(initial_penetration*cos(constants.beta*current_theta) +...
        (initial_theta/constants.beta + initial_penetration) *sin(constants.beta*current_theta));...
        current_ind_tyre = current_ind_tyre + direction;
    current_theta = current_theta + dTheta;
end
end
% State transition functions, functions with in-place write that change
% objects states

function road = createRoad(road, constants)
step_x = linspace(0 , constants.step_width, 100)' ;
step_y = constants.step_height * sin(step_x *constants.step_phase_width/constants.step_width - pi/2);
step_start_ind = floor((constants.step_position - constants.step_width/2) / constants.road_distance_step);
step_end_ind = ceil(constants.step_width/constants.road_distance_step) + step_start_ind;
step_y = road.y(step_end_ind) + step_y - step_y(1);

road.x = [road.x(1:step_start_ind-1);step_x + road.x(step_start_ind); road.x(step_end_ind+1:end)];
road.y = [road.y(1:step_start_ind-1) ; step_y; road.y(step_end_ind +1:end) - road.y(step_end_ind+1) + step_y(end)];

road.gradient = getCentralDiffDerivative(road.x , road.y);
road.ddy = getCentralDiffDerivative(road.x , road.gradient);
road.curvature = road.ddy./((1 + road.gradient.^2).^(3/2));
end

function [tyre, road] = getAllPenetrations(tyre, road)
search_start_ind = find(road.x < tyre.x_centre - tyre.free_radius , 1 , "last");
search_end_ind = find(road.x > tyre.x_centre +tyre.free_radius, 1 , "first");
if(isempty(search_start_ind))
    search_start_ind = 1;
end
if(isempty(search_end_ind))
    search_end_ind = 1;
end
road.penetration_inds = [];
new_penetration_inds = findPenetratingInds(tyre.x_centre, tyre.y_centre, road.x, road.y,...
    tyre.free_radius, search_start_ind, search_end_ind);
while(~isempty(new_penetration_inds))
    road.penetration_inds(end+1 , :) = new_penetration_inds;
    new_penetration_inds = findPenetratingInds(tyre.x_centre, tyre.y_centre, road.x, road.y,...
        tyre.free_radius, new_penetration_inds(2), search_end_ind);
end
tyre = getTyrePenetrationInds(tyre , road);
end

function [tyre, road] = getContactCentres(tyre ,road)
road.contact_centre_inds = zeros(length(road.penetration_inds(: , 1)), 1);
tyre.contact_centre_inds  = zeros(length(road.penetration_inds(: , 1)),1);
for i = 1: length(tyre.penetration_inds(: , 1))
    x = road.x(road.penetration_inds(i , 1) : road.penetration_inds(i , 2)) - tyre.x_centre;
    y = road.y(road.penetration_inds(i , 1) : road.penetration_inds(i , 2)) - tyre.y_centre;
    distances = sqrt(x.^2 + y.^2);
    [~ , min_ind] = min(distances);
    road.contact_centre_inds(i) = road.penetration_inds(i , 1) + min_ind -1;
    tyre.contact_centre_inds(i) = getTyreIndexfromRoadIndex(tyre , road, road.contact_centre_inds(i));
end
end

function tyre = getTyrePenetrationInds(tyre, road)
tyre.penetration_inds = zeros(size(road.penetration_inds));
for i = 1:length(tyre.penetration_inds( : , 1))
    tyre.penetration_inds(i  ,1) = getTyreIndexfromRoadIndex(tyre, road,  road.penetration_inds(i , 1));
    tyre.penetration_inds(i  ,2) = getTyreIndexfromRoadIndex(tyre, road, road.penetration_inds(i , 2));
end
end

function [tyre, road] = getContactBoundaries(tyre, road, constants)
% TODO, implement the actual separation condition here later, for now
% it only uses the different between derivatives.
for i = 1:length(road.contact_centre_inds)
    % find aft boundary
    for check_ind = road.contact_centre_inds(i):-1:road.penetration_inds(i , 1)
        road.contact_boundary_inds(i , 1) = check_ind;
        tyre.contact_boundary_inds(i, 1) = getTyreIndexfromRoadIndex(tyre, road , check_ind);
        %         separation_angle = pi - atan2(road.y(check_ind-1) - road.y(road.contact_centre_inds) , ...
        %             road.x(check_ind-1) - road.x(road.contact_centre_inds));
        separation_angle = road.gradient(check_ind) - road.gradient(road.contact_centre_inds(i))
        if separation_angle > constants.separation_threshold
            break;
        end
    end
    % find fore boundary
    for check_ind = road.contact_centre_inds(i):1:road.penetration_inds(i , 2)
        road.contact_boundary_inds(i , 2) = check_ind;
        tyre.contact_boundary_inds(i, 2) = getTyreIndexfromRoadIndex(tyre, road , check_ind);
        separation_angle = atan2(road.y(check_ind) - road.y(road.contact_centre_inds) , ...
            road.x(check_ind) - road.x(road.contact_centre_inds));
        if separation_angle > constants.separation_threshold
            break;
        end
    end
end
end

% plotting functions

function plot_handles = drawRoad(road, plot_handles, fig_handle)
if nargin == 2
    fig_handle = gcf();
end
figure(fig_handle);
hold on
plot_handles(end+1) = plot(road.x, road.y , 'black');
end

function plot_handles = drawUndeformedTyre(tyre,plot_handles, fig_handle)
if nargin == 2
    fig_handle = gcf();
end
x  = tyre.free_radius*cos(tyre.theta) + tyre.x_centre;
y = tyre.free_radius* sin(tyre.theta) + tyre.y_centre;
figure(fig_handle);
x = [x;x(1)];
y = [y;(1)];
plot_handles(end+1) = plot(x , y , 'r.-', MarkerSize=10);
end

function plot_handles = drawPenetrations(tyre, road , plot_handles ,fig_handle)
if nargin == 3
    fig_handle = gcf();
end
figure(fig_handle);
hold on
for i = 1:length(tyre.penetration_inds(: , 1))
    plot_handles(end+1) = plot(road.x(road.penetration_inds(i , 1) : road.penetration_inds(i, 2)),...
        road.y(road.penetration_inds(i , 1):road.penetration_inds(i , 2)),".", MarkerSize=5);
    plot_handles(end+1) =plot(road.x(road.contact_centre_inds(i)), road.y(road.contact_centre_inds(i)), "ro", MarkerSize=10);
    plot_handles(end+1) = plot(tyre.x_centre + tyre.free_radius*cos(tyre.theta(tyre.penetration_inds(i , 1))), ...
        tyre.y_centre + tyre.free_radius*sin(tyre.theta(tyre.penetration_inds(i , 1))), "o", MarkerSize=10);
    plot_handles(end+1) = plot(tyre.x_centre + tyre.free_radius*cos(tyre.theta(tyre.penetration_inds(i , 2))), ...
        tyre.y_centre + tyre.free_radius*sin(tyre.theta(tyre.penetration_inds(i , 2))), "o", MarkerSize=10);
    plot_handles(end+1) = plot(tyre.x_centre + tyre.free_radius * cos(tyre.theta(tyre.contact_centre_inds)), ...
        tyre.y_centre + tyre.free_radius * sin(tyre.theta(tyre.contact_centre_inds)), "m.", MarkerSize= 10);
    plot_handles(end+1) = plot(road.x(road.contact_boundary_inds(i , 1)), ...
        road.y(road.contact_boundary_inds(i , 1)), "g*", MarkerSize= 15);
    plot_handles(end+1) = plot(road.x(road.contact_boundary_inds(i , 2)),...
        road.y(road.contact_boundary_inds(i , 2)), "g*", MarkerSize= 15);
end
end

