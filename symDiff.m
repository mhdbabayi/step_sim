cos(alpha)*((D*cos(theta))/(R*(1 - (D^2*sin(theta)^2)/R^2)^(1/2)) - 1)

%%
clear
clc
syms D R r theta

alpha = asin((D/R)*sin(theta)) - theta;
w0 = r - R*(sin(alpha)/sin(theta));
d1 = diff(sin(alpha) , theta);
subs(d1 , asin(D*sin(theta)/R) - theta, 'alpha')
%%
syms theta tyre_radius R tyre_centre_height  beta theta0 gamma r G
syms theta 
% tyre_radius = 0.33;
% terrain_radius = 0.2;
% tyre_centre_height = tyre_radius*(1 - 0.6);
% beta = 5;
u = (1 + tyre_centre_height/R);
alpha = asin(u*sin(theta)) - theta;
max_theta = acos(((tyre_centre_height+R)^2 + tyre_radius^2 - R^2)/(2*(tyre_centre_height+R)*tyre_radius));
max_theta_func = matlabFunction(max_theta);
theta_range = linspace(0 , max_theta, 100);
r = sqrt((tyre_centre_height+ R)^2 + R^2 - 2*(tyre_centre_height+R)*R*cos(alpha));
dr_dtheta = diff(alpha , theta)*R*sin(alpha + theta);
ddr_dtheta = diff(dr_dtheta , theta);
dr0 = tyre_radius - r;
dr_dtheta0 = -dr_dtheta;
gamma = theta - theta0;
beam = exp(-beta* (gamma))*((dr0 - dr_dtheta0/beta)*sin(beta*(gamma)) + (-dr0)*cos(beta*(gamma)));
d_beam_dtheta = diff(beam , theta);
dd_beam_dtheta = diff(d_beam_dtheta, theta);
dd_beam_func = matlabFunction(subs(dd_beam_dtheta , theta0 , theta));
dd_road_func = matlabFunction(-ddr_dtheta);
tyre_radius = 0.77/2;
beta = 10;
terrain_radius_range = [0.01:0.001:0.2 ];
penetration_range = [0.005 :0.01:0.15];

results = table([], [], []);
results.Properties.VariableNames = {'terrain_radius', 'penetration' , 'separation'};
for terrain_radius_num = terrain_radius_range
    for tyre_centre_height_num = tyre_radius - penetration_range
        theta_numerical = [0:deg2rad(0.1) , max_theta_func(...
            terrain_radius_num,tyre_radius,tyre_centre_height_num)];
%         eval(subs(theta_range ,...
%             {terrain_radius , tyre_centre_height , tyre_radius},...
%             {terrain_radius_num, tyre_centre_height_num, tyre_radius} ));
        lhs = dd_beam_func(beta , terrain_radius_num, theta_numerical, tyre_radius, tyre_centre_height_num);
        rhs = dd_road_func(terrain_radius_num,theta_numerical, tyre_centre_height_num);
        crossing = theta_numerical(find(lhs > rhs , 1 , 'first'));
        if isempty(crossing)
            crossing = theta_numerical(end);
        end
        results(end + 1 , :) = {terrain_radius_num ,tyre_radius - tyre_centre_height_num , crossing};
        plot(rad2deg(theta_numerical) , lhs);
        hold on 
        plot(rad2deg(theta_numerical) , rhs , "*");
    end
    fprintf("terrain radius %.1fmm done\n", 1000*terrain_radius_num);
end
ylim([-10 , 10])

%%
close all
subplot(1 , 2 , 1);
clc
unique_penetration = unique(results.penetration);
legend_labels = {};
for h = 1:length(unique_penetration)
    h_idx = find(results.penetration == unique_penetration(h));
    y = results.separation(h_idx) * tyre_radius*1000;
    x = results.terrain_radius(h_idx)*1000;
    plot(x , y);
    plot_text = "pen: " + num2str(unique_penetration(h)*1000)+ "mm";
    text_pos_idx = length(x);
    if ~isempty(text_pos_idx)
        text(x(text_pos_idx) , y(text_pos_idx) , plot_text);
    end
    hold on
end
xlabel("terrain radius");
ylabel("separation mm");
% daspect([1 , 1, 1])
% subplot(1 , 2, 2)
% unique_terrain_radius = unique(results.terrain_radius);
% legend_labels = {};
% for h = 1:length(unique_terrain_radius)
%     h_idx = find(results.tyre_height == unique_terrain_radius(h));
%     plot(results.(h_idx) , rad2deg(results.separation(h_idx)));
%     legend_labels{end+1} = "penetration: " + num2str(tyre_radius -unique_penetration(h));
%     hold on
% end
% legend(legend_labels);
% xlabel("terrain radius");
% ylabel("separation");
%%
lhs = eval(subs(dd_beam_dtheta,{gamma,theta}, {0,theta_range}));
rhs = eval(subs(-ddr_dtheta , theta , theta_range));
plot(rad2deg(theta_range) , lhs , '*');
hold on 
plot(rad2deg(theta_range) , rhs);
[tol , solution_idx ]= min(abs(lhs - rhs));
solution_theta = theta_range(solution_idx);
plot(rad2deg(solution_theta) , lhs(solution_idx), "m*")
theta_range = linspace(solution_theta - deg2rad(1) , solution_theta+deg2rad(1), 100);
lhs = eval(subs(dd_beam_dtheta,{gamma,theta}, {0,theta_range}));
rhs = eval(subs(-ddr_dtheta , theta , theta_range));
plot(rad2deg(theta_range) , lhs , '*');
hold on 
plot(rad2deg(theta_range) , rhs);
%%
% clear
clc
close all
tyre_radius = 0.33;
R = 1;
tyre_centre_height = tyre_radius*(1 - 0.1);
beta = 5;
sep_theta = [];
penetration = linspace(0 , 0.1, 100);
for tyre_centre_height = (tyre_radius - penetration)
sep_theta(end+1) = get_separation(beta , tyre_radius , R , tyre_centre_height);
end
figure
plot(penetration , rad2deg(sep_theta(: , 1)));
xlabel ("tyre penetration")
ylabel ("separation angle")
%%
function separation_theta = get_separation(beta,tyre_radius, terrain_radius, tyre_centre_height)
     syms theta theta0
     u = (1 + tyre_centre_height/terrain_radius);
    alpha = asin(u*sin(theta)) - theta;
    max_theta = acos(((tyre_centre_height+terrain_radius)^2 + tyre_radius^2 - terrain_radius^2)/(2*(tyre_centre_height+terrain_radius)*tyre_radius));
    theta_range = [0:deg2rad(0.5):max_theta];
    r = sqrt((tyre_centre_height+ terrain_radius)^2 + terrain_radius^2 - 2*(tyre_centre_height+terrain_radius)*terrain_radius*cos(alpha));
    dr_dtheta = diff(alpha , theta)*terrain_radius*sin(alpha + theta);
    ddr_dtheta = diff(dr_dtheta , theta);
    dr0 = tyre_radius - r;
    dr_dtheta0 = -dr_dtheta;
    gamma = theta - theta0;
    beam = exp(-beta* (gamma))*((dr0 - dr_dtheta0/beta)*sin(beta*(gamma)) + (-dr0)*cos(beta*(gamma)));
    d_beam_dtheta = diff(beam , theta);
    dd_beam_dtheta = diff(d_beam_dtheta, theta);  
    lhs = eval(subs(dd_beam_dtheta,{gamma,theta}, {0,theta_range}));
    rhs = eval(subs(-ddr_dtheta , theta , theta_range));
    plot(rad2deg(theta_range) , lhs , '*');
    hold on 
    plot(rad2deg(theta_range) , rhs);
    solution_idx= find(lhs > rhs , 1 , 'first');
    if isempty(solution_idx)
        solution_idx = length(lhs);
    end

    solution_theta = theta_range(solution_idx);
    plot(rad2deg(solution_theta) , lhs(solution_idx), "m.", MarkerSize=30)
%     theta_range = [solution_theta -deg2rad(1):deg2rad(0.1):solution_theta + deg2rad(1)];
%     lhs = eval(subs(dd_beam_dtheta,{gamma,theta}, {0,theta_range}));
%     rhs = eval(subs(-ddr_dtheta , theta , theta_range));
%     plot(rad2deg(theta_range) , lhs , '.');
%     hold on 
%     plot(rad2deg(theta_range) , rhs);
    separation_theta = solution_theta;
end
