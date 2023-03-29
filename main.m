clear
clc
close all
tyre_sim_structs
road_s = createRoad(road_s, constants_s);
main_fig= figure();
xlim([road_s.x(1) , road_s.x(end)]);
ylim([tyre_s.y_centre - tyre_s.free_radius*1.1 , tyre_s.y_centre + tyre_s.free_radius*1.2]);
hold on
daspect([1 1 1])
plot_handles_s = drawRoad(road_s, plot_handles);
% vidRecorder = VideoWriter("beamVid");
% vidRecorder.open();
%%
clear road tyre constants
while tyre_s.x_centre < (constants_s.step_position + constants_s.tyre_radius)
    plot_handles_s = drawUndeformedTyre(tyre_s, plot_handles_s);
    try
    [tyre_s, road_s] = tyreRoadContact(tyre_s , road_s , constants_s);
    catch me
        disp(me)
        tyre_s.x_centre = tyre_s.x_centre + 0.01;
        continue;
    end
    plot_handles_s = drawPenetrations(tyre_s, road_s , plot_handles_s);
    tyre_s = getDeformedProfile(tyre_s ,road_s ,constants_s , length(tyre_s.contact_centre_inds));
    plot_handles_s = drawDeformedTyre(tyre_s , plot_handles_s);
    pause();
    delete(plot_handles_s.frame);
    plot_handles_s.frame = [];
    tyre_s.x_centre = tyre_s.x_centre + 0.01;
end
% vidRecorder.close();

%% testing second derivative
% inds = (road_s.contact_centre_inds(1):road_s.penetration_inds(1 , 2))';
% ddr = polarSecondDerivative(road_s.x(inds) - tyre_s.x_centre , road_s.y(inds) - tyre_s.y_centre , road_s.gradient(inds), road_s.ddy(inds));
% [ddr , road2tyreDerivativeTF()]
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

function road_inds = getRoadIndexfromTyreIndex(tyre , road , tyre_ind)
    x = tyre.x_centre + tyre.r_free(tyre_ind).*cos(tyre.theta(tyre_ind));
    road_inds = interp1(road.x , (1:length(road.x)), x , "nearest");
end

function polar_coordinates = road2tyreTF(tyre, road_point_x, road_point_y)
polar_coordinates(1 ,1) = sqrt((road_point_x - tyre.x_centre)^2 + (road_point_y - tyre.y_centre)^2);
polar_coordinates(2 ,1) = mod(atan2(road_point_y - tyre.y_centre , road_point_x - tyre.x_centre), 2*pi);
end

function polar_derivative = road2tyreDerivativeTF(tyre, road, road_point_ind)
x = road.x(road_point_ind) -tyre.x_centre;
y = road.y(road_point_ind) - tyre.y_centre;
R = sqrt(x^2 + y^2);
theta = atan2(y , x);
y_prime = road.gradient(road_point_ind);
polar_derivative = (x + y*y_prime)/(y_prime * cos(theta) - sin(theta));
end

function polar_derivative = polarDerivative(x , y , dy)
R = sqrt(x.^2 + y.^2);
polar_derivative = R.*(x + y.*dy)./(x.*dy - y);
end

function polar_second_derivative  = polarSecondDerivative(x , y , dy , ddy)
r = sqrt(x.^2 + y.^2);
theta = atan2(y , x);
dr_dtheta = polarDerivative(x , y , dy);
dx_dtheta = dr_dtheta.*x./r  - y;
dy_dtheta = dr_dtheta.*y./r + x;
% dydx_dtheta = ddy .* dx_dtheta + (ddy./dy).*dy_dtheta; % d(dy/dx)/dtheta;
dydx_dtheta = ddy .*dx_dtheta;
polar_second_derivative = ((x.*dy - y).*(dx_dtheta + y.*dydx_dtheta + dy.*dy_dtheta) - ...
    (y.*dy + x).*(dy.*dx_dtheta + dydx_dtheta.*x - dy_dtheta)) .* r./ ...
    (x.*dy - y).^2 +  dr_dtheta.*(y.*dy  +x)./(x.*dy - y);

end

function global_coordinates = tyre2roadTF(tyre, theta , r)
global_coordinates(1 , 1) = tyre.x_centre + r*cos(theta);
global_coordinates(2 , 1) = tyre.y_centre + r*sin(theta);
end

function initial_conditions = getInitialConditions(tyre, road , contact_ind)
initial_conditions.fore = [];
initial_conditions.aft = [];
road_ind_fore = road.contact_boundary_inds(contact_ind , 2);
road_ind_aft = road.contact_boundary_inds(contact_ind , 1);
fore_point_polar = road2tyreTF(tyre , road.x(road_ind_fore),road.y(road_ind_fore));
aft_point_polar =  road2tyreTF(tyre , road.x(road_ind_aft),road.y(road_ind_aft));
initial_conditions.fore.w_prime = -road2tyreDerivativeTF(tyre ,road , road_ind_fore);
initial_conditions.aft.w_prime = road2tyreDerivativeTF(tyre , road, road_ind_aft);
initial_conditions.fore.w =   tyre.free_radius - fore_point_polar(1);
initial_conditions.aft.w =   tyre.free_radius - aft_point_polar(1);
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
% road_contact_inds = (road.contact_centre_inds(1):road.penetration_inds(1 , 2))';
% ddr_dtheta = polarSecondDerivative(road_s.x(road_contact_inds) - tyre_s.x_centre , ...
%     road_s.y(road_contact_inds) - tyre_s.y_centre ,...
%     road_s.gradient(road_contact_inds), road_s.ddy(road_contact_inds));
% dr_dtheta = polarDerivative(road_s.x(road_contact_inds) - tyre_s.x_centre , ...
%     road_s.y(road_contact_inds) - tyre_s.y_centre ,...
%     road_s.gradient(road_contact_inds));
% road.contact_boundary_inds(1 , 2) = find(2 * constants.beta*dr_dtheta > ddr_dtheta , 1 , 'first') + road.contact_centre_inds(1);
% if isempty( road.contact_boundary_inds(1 ,2))
%     road.contact_boundary_inds(1 , 2) = road.penetration_inds(1 ,2);
% end
% 
% road_contact_inds = (road.contact_boundary_inds(1,1):road.contact_centre_inds(1))';
% ddr_dtheta = polarSecondDerivative(road_s.x(road_contact_inds) - tyre_s.x_centre , ...
%     road_s.y(road_contact_inds) - tyre_s.y_centre ,...
%     road_s.gradient(road_contact_inds), road_s.ddy(road_contact_inds));
% dr_dtheta = polarDerivative(road_s.x(road_contact_inds) - tyre_s.x_centre , ...
%     road_s.y(road_contact_inds) - tyre_s.y_centre ,...
%     road_s.gradient(road_contact_inds));
% road.contact_boundary_inds(1 , 2) = find(2 * constants.beta*dr_dtheta > ddr_dtheta , 1 , 'first') + road.contact_centre_inds(1);
% if isempty( road.contact_boundary_inds(1 ,2))
%     road.contact_boundary_inds(1 , 2) = road.penetration_inds(1 ,2);
% end

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

function tyre = getDeformedProfile(tyre, road, constants , contact_ind)
initial_conditions = getInitialConditions(tyre, road, contact_ind);
dTheta = 2*pi/(length(tyre.theta));
theta = (0:dTheta:pi/2)';
% fore
deflection_profile = exp(-constants.beta * theta).*(initial_conditions.fore.w*cos(constants.beta*theta) +...
    (initial_conditions.fore.w_prime/constants.beta + initial_conditions.fore.w) *sin(constants.beta*theta));
shifted_inds = mod((0:(length(deflection_profile)-1))'+ tyre.contact_boundary_inds(contact_ind , 2) -1 , length(tyre.theta)) + 1;
tyre.r_deformed(shifted_inds) = tyre.r_free(shifted_inds) - deflection_profile;
%aft
deflection_profile = exp(-constants.beta * theta).*(initial_conditions.aft.w*cos(constants.beta*theta) +...
    (initial_conditions.aft.w_prime/constants.beta + initial_conditions.aft.w) *sin(constants.beta*theta));
shifted_inds = mod(-(0:(length(deflection_profile)-1))' + tyre.contact_boundary_inds(contact_ind , 1) -1 , length(tyre.theta)) + 1;
tyre.r_deformed(shifted_inds) = tyre.r_free(shifted_inds) - deflection_profile;
% TODO 
% make it go round
tyre_contact_inds = (tyre.contact_boundary_inds(contact_ind , 1):tyre.contact_boundary_inds(contact_ind , 2))';
road_contact_inds = getRoadIndexfromTyreIndex(tyre , road, tyre_contact_inds);
tyre.r_deformed(tyre_contact_inds) = sqrt((road.x(road_contact_inds) - tyre.x_centre).^2 + (road.y(road_contact_inds) - tyre.y_centre).^2);
end

% plotting functions

function plot_handles = drawRoad(road, plot_handles, fig_handle)
if nargin == 2
    fig_handle = gcf();
end
figure(fig_handle);
hold on
plot_handles.static(end+1) = plot(road.x, road.y , 'black');
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
plot_handles.frame(end+1) = plot(x , y , 'r.-', MarkerSize=3);
end

function plot_handles = drawDeformedTyre(tyre,plot_handles, fig_handle)
if nargin == 2
    fig_handle = gcf();
end
x  = tyre.r_deformed.*cos(tyre.theta) + tyre.x_centre;
y = tyre.r_deformed.* sin(tyre.theta) + tyre.y_centre;
figure(fig_handle);
x = [x;x(1)];
y = [y;(1)];
plot_handles.frame(end+1) = plot(x , y , 'm.-', MarkerSize=5);
end

function plot_handles = drawPenetrations(tyre, road , plot_handles ,fig_handle)
if nargin == 3
    fig_handle = gcf();
end
figure(fig_handle);
hold on
for i = 1:length(tyre.penetration_inds(: , 1))
    plot_handles.frame(end+1) = plot(road.x(road.penetration_inds(i , 1) : road.penetration_inds(i, 2)),...
        road.y(road.penetration_inds(i , 1):road.penetration_inds(i , 2)),".", MarkerSize=5);
    plot_handles.frame(end+1) =plot(road.x(road.contact_centre_inds(i)), road.y(road.contact_centre_inds(i)), "ro", MarkerSize=10);
    plot_handles.frame(end+1) = plot(tyre.x_centre + tyre.free_radius*cos(tyre.theta(tyre.penetration_inds(i , 1))), ...
        tyre.y_centre + tyre.free_radius*sin(tyre.theta(tyre.penetration_inds(i , 1))), "o", MarkerSize=10);
    plot_handles.frame(end+1) = plot(tyre.x_centre + tyre.free_radius*cos(tyre.theta(tyre.penetration_inds(i , 2))), ...
        tyre.y_centre + tyre.free_radius*sin(tyre.theta(tyre.penetration_inds(i , 2))), "o", MarkerSize=10);
    plot_handles.frame(end+1) = plot(tyre.x_centre + tyre.free_radius * cos(tyre.theta(tyre.contact_centre_inds)), ...
        tyre.y_centre + tyre.free_radius * sin(tyre.theta(tyre.contact_centre_inds)), "m.", MarkerSize= 10);
    plot_handles.frame(end+1) = plot(road.x(road.contact_boundary_inds(i , 1)), ...
        road.y(road.contact_boundary_inds(i , 1)), "g*", MarkerSize= 15);
    plot_handles.frame(end+1) = plot(road.x(road.contact_boundary_inds(i , 2)),...
        road.y(road.contact_boundary_inds(i , 2)), "g*", MarkerSize= 15);
end
end

