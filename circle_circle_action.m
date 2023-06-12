%%
clear
close all
clc
%%%Terrain is a circle with radius terrain_radius, and tyre is a circle with radius
%%%tyre_radius and distance D from from the centre of terrain
%%% r (lower_case) refers to a vector from tyre centre to the surface of
%%% the terrain
% %% R is a vector from a point on the surface of the terrain to the centre of the terrain
% %% w is the amount of deflection at a point on the
    %%% circumference of the deflected tyre (positive in compression)
%%% beta is the a function of bending stiffness, second moment of area of
%%% cross section and compression stiffness in the tyre, explained more
%%% clearly in the beam solution.
% %% for any point in space:
%%%         theta is angle between a vector from tyre centre and that point
%%%         and concentric line, measure from right
%%%         alpha is angle between a vector from terrain centre and that
%%%         point, relative to concentric line measured from left
%%% we aim to find the point where (d^2/dtheta^2) of w, assuming that point
%%% is the initial point of the beam solution, is greater than
%%% (d^2/dtheta^2) of R
%%% anything with _v at the end is a vector with origin at tyre centre
clear
clc
close all
syms theta theta0 R penetration w r terrain_radius beta
tyre_radius = 0.788/2;
D = terrain_radius + tyre_radius - penetration;
tyre_height = tyre_radius - penetration;
alpha = asin((1+tyre_height/terrain_radius)*sin(theta0)) - theta0;
terrain_centre_v = [D;0];
R_v = terrain_centre_v + terrain_radius*[-cos(alpha);sin(alpha)];
R = sqrt(dot(R_v , R_v));
dR_dtheta = diff(tyre_radius -R , theta0);
ddR_dtheta_func = matlabFunction(real(diff(dR_dtheta , theta0)));
%%% delta_w at the initial_conditions
w0 = tyre_radius - R;
dw0 = diff(w0, theta0);
w = exp(-beta * (theta - theta0))*(...
    (w0 + dw0/beta)*sin(beta*(theta- theta0)) + ...
    w0 * cos(beta*(theta - theta0)));
dw_dtheta = diff(w , theta);
ddw_dtheta_func = matlabFunction(real(diff(dw_dtheta , theta)));
%%% law of cosines
max_theta_func = matlabFunction(acos((tyre_radius^2 + D^2 - terrain_radius^2)/(2*D*tyre_radius)));
beta0 = 5;
%% testing 
close all
clc
pen = 0.009;
tr = 0.5;
beta0 = 5;
th = get_theta_range(tr , tyre_radius , pen, 0.01);
ddw = ddw_dtheta_func(beta0, pen , tr ,th, th );
ddr = ddR_dtheta_func(pen, tr , th);
plot(rad2deg(th) , ddw, "--");
hold on 
plot(rad2deg(th) , ddr);
figure
sep_idx = get_separation_idx(ddw , ddr , th);
plot(rad2deg(th(sep_idx)), ddw(sep_idx), "*", MarkerSize=15);
terrain_points_angles = linspace(0 , 2*pi ,3600);
terrian_points = 1000 * tr*[cos(terrain_points_angles);sin(terrain_points_angles)];

tyre_point_angles = terrain_points_angles + 3*pi/2 + th(sep_idx);
bc1 = eval(subs(w0, symvar(w0) , {pen , tr , th(sep_idx)}));
bc2 = bc1 + eval(subs(dw0, symvar(dw0) , {pen , tr , th(sep_idx)}))/beta0;
w_angles = tyre_point_angles - tyre_point_angles(1);
tyre_w = tyre_radius - exp(-beta0 * (w_angles)).*(...
    bc1* cos(beta0 * w_angles) + bc2*sin(beta0*w_angles));
tyre_points = 1000*tyre_w .* [cos(tyre_point_angles); sin(tyre_point_angles)] + 1000*[0;tr + tyre_radius - pen];
hold on 
plot(tyre_points(1 , :) , tyre_points(2 , :) , "*-")
plot(terrian_points(1 , :) , terrian_points(2 , :))
daspect([1 , 1 ,1])
legend({"beam", "terrain"})
%%
clc
beta0 = 5;
results = table([], [], []);
results.Properties.VariableNames = {'terrain_radius', 'penetration' , 'separation'};
penetration_num = [0.001:0.002:0.25];
terrain_radius_num = [0.01:0.002:2];
results = zeros(length(penetration_num),length(terrain_radius_num));
for pen_idx = 1:length(penetration_num)
    for terrain_idx = 1:length(terrain_radius_num)
        this_penetration = penetration_num(pen_idx);
        this_terrain_radius = terrain_radius_num(terrain_idx);
        if (this_penetration < 0.02) 
            angle_step_deg = 0.01;
        else
            angle_step_deg = 1;
        end
        theta_num =get_theta_range(this_terrain_radius, tyre_radius , this_penetration, angle_step_deg);
        ddw_num = ddw_dtheta_func(beta0,this_penetration,...
                                                        this_terrain_radius,theta_num , theta_num);
        ddR_num = ddR_dtheta_func(this_penetration, this_terrain_radius, theta_num);
        sep_idx = get_separation_idx(ddw_num , ddR_num , theta_num);
        results(pen_idx , terrain_idx) = theta_num(sep_idx);
    end
    fprintf("penetration %.1fmm done\n", penetration_num(pen_idx)*1000)
end
%% plotting all penetrations for a single terrain radius
close all
terrain_R_mm = 500;
for terrain_R_mm = [10:20:1000]
tr_idx = find(terrain_radius_num*1000 > terrain_R_mm , 1 , 'first');
y = tyre_radius*1000*results(: , tr_idx);

plot(penetration_num*1000, y)
hold on
text(penetration_num(end)*1000,y(end), sprintf("R: %.0fmm", terrain_radius_num(tr_idx)*1000))
end
xlabel("penetration mm")
ylabel("separation theta deg\circ")

daspect([1 , 1 ,1])
%% get fitted objects for each penetration
clc
close all
fitted_results = results;
fit_coefficients = nan(length(penetration_num) , 4);
for pen_idx = 1:length(penetration_num)
y  = results(pen_idx , :);
x = terrain_radius_num * 1000;
start_idx = find(y > 0 , 1 , 'first');
if ~isempty(start_idx)
fit_x = x(start_idx:end);
fit_y = y(start_idx:end);
fit_opts = fitoptions('rat11', 'StartPoint' , [fit_y(end) , 0 , 0]);
fit_object = fit(fit_x' , fit_y' , 'rat11', fit_opts);
fit_coefficients(pen_idx ,:) = [coeffvalues(fit_object), fit_x(start_idx)];
fitted_results(pen_idx , start_idx:end) = fit_object(fit_x);
end
end
clc
close all
figure
surf(terrain_radius_num*1000 , penetration_num*1000, tyre_radius*1000*results, EdgeColor="none")
hold on 
surf(terrain_radius_num*1000 , penetration_num*1000, tyre_radius*1000*fitted_results, EdgeColor="r")
xlabel("terrain radius mm");
ylabel("penetration mm")
zlabel("separation angle deg\circ")
%%
close all
clc
penetration_mm = 9;
pen_idx = find(penetration_num*1000 > penetration_mm, 1 , 'first');
batch_fitted_y  = fitted_results(pen_idx , :);
original_y = results(pen_idx , :);
x = terrain_radius_num * 1000;
% fit_obj = fit(x(original_y>0)' , original_y(original_y>0)' , 'rat11');

plot(x, 1000*tyre_radius*batch_fitted_y, "b*-");
hold on 
plot(x , tyre_radius*1000*original_y, LineWidth=5)
% plot(x(original_y>0) , 1000*tyre_radius*fit_obj(x(original_y> 0)) , "mx")

xlabel("terrain radius mm")
ylabel("contact length mm")
text(x(end), 1000*tyre_radius*original_y(end), sprintf("pen: %.0fmm", penetration_num(pen_idx)*1000))
hold on

% daspect([ 1 , 1 ,1])
%%

%%

%% helper functions 
function theta_range = get_theta_range(terrain_radius, tyre_radius , penetration, angle_step_deg)
centre_distance = tyre_radius + terrain_radius - penetration;    
if penetration > terrain_radius
    max_theta = atan(terrain_radius/centre_distance);
else
    max_theta = acos((tyre_radius^2 + centre_distance^2 - terrain_radius^2)...
        /(2*tyre_radius*centre_distance));
end
     theta_range = (0:deg2rad(angle_step_deg): max_theta);
end


function separation_idx = get_separation_idx(ddw, ddr, theta)
separation_idx = find(ddw > ddr , 1 , 'first');
if isempty(separation_idx)
    separation_idx = length(theta);
end
end
