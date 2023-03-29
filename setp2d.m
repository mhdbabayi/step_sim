clear
clc
close all
%% constant parameters
tyre_radius = 1;
tyre_max_enveloping_angle = deg2rad(20);
total_length = tyre_radius * 5;
step_position = tyre_radius * 3;
num_step_points = 100;
path_resolution = 0.01;
step_height = tyre_radius/3;
step_width = tyre_radius*0.1;
%% Initializing objects
%step
step_y = step_height*(sin(linspace(-pi/2, pi/2, num_step_points)')+ 1)/2;
step_x = linspace(step_position - step_width/2, step_position + step_width/2 , num_step_points)';
%path
path.x = [(0:path_resolution:step_x(1) - path_resolution)';...
    step_x;(step_x(end)+path_resolution:path_resolution:total_length)'];
path.y= [(0:path_resolution:step_x(1) - path_resolution)'*0 + step_y(1);step_y;...
    (step_x(end)+path_resolution:path_resolution:total_length)'*0 + step_y(end)];
path.tangent = atan([0;diff(path.y)./diff(path.x)]);
ds = sqrt(diff(path.y).^2 + diff(path.x).^2);
path.ds = [ds;ds(end)];
path.curvature = [0;diff(path.tangent)]./path.ds;
%tyre
tyre_position.x = tyre_radius*2.3;
tyre_position.y = tyre_radius*0.9 + path.y(1);
%% plotting
close all;
fig = figure();
plot(path.x , path.y , '.-');
hold on
viscircles([tyre_position.x , tyre_position.y], tyre_radius);
xlim([0 , total_length]);
ylim([-tyre_radius/2 ,1.5*tyre_radius])
daspect([1,1,1])
%% finding contacts
penetrations = getPenetration(tyre_position.x, tyre_position.y , tyre_radius, path.x, path.y);
for pen = penetrations
    plot(path.x([pen.start, pen.stop]), path.y([pen.start , pen.stop]), '*');
    plot(path.x(pen.closest_ind) , path.y(pen.closest_ind),  'o');
    [envelop_start,envelop_stop] = getContact(path, pen.start, pen.stop, pen.closest_ind, deg2rad(20));

    [circle_centre, circle_radius, start_angle, end_angle]  = getDeformedSection(tyre_position.x , tyre_position.y,...
        tyre_radius, path.x(envelop_stop), path.y(envelop_stop), path.tangent(envelop_stop+1));
    %viscircles(circle_centre' ,abs(circle_radius) , 'Color','bla');

    start_point = circle_centre + circle_radius*([cos(start_angle); sin(start_angle)]);
    %plot([circle_centre(1), start_point(1)] , [circle_centre(2), start_point(2)], 'm')
    end_point = circle_centre + circle_radius*([cos(end_angle); sin(end_angle)]);
    %plot([circle_centre(1), end_point(1)] , [circle_centre(2), end_point(2)], 'm')

    circle_angles  = linspace(min(start_angle , end_angle), max(start_angle , end_angle), 100);
    circle_points = [[circle_centre(1) + circle_radius*cos(circle_angles)]',...
        [circle_centre(2) + circle_radius*sin(circle_angles)]'];
    plot(circle_points(: , 1) , circle_points(: , 2) ,'bla', LineWidth=3)
    [circle_centre, circle_radius, start_angle, end_angle] = getDeformedSection(tyre_position.x , tyre_position.y,...
        tyre_radius, path.x(envelop_start), path.y(envelop_start), path.tangent(envelop_start-1));
    %viscircles(circle_centre' ,abs(circle_radius), 'Color','b');
    start_point = circle_centre + circle_radius*([cos(start_angle); sin(start_angle)]);
    %plot([circle_centre(1), start_point(1)] , [circle_centre(2), start_point(2)], 'g')
    end_point = circle_centre + circle_radius*([cos(end_angle); sin(end_angle)]);
    %plot([circle_centre(1), end_point(1)] , [circle_centre(2), end_point(2)], 'g')

    circle_angles  = linspace(min(start_angle , end_angle), max(start_angle , end_angle), 100);
    circle_points = [[circle_centre(1) + circle_radius*cos(circle_angles)]',...
        [circle_centre(2) + circle_radius*sin(circle_angles)]'];
    plot(circle_points(: , 1) , circle_points(: , 2) ,'bla' , LineWidth=3)

    plot(path.x(envelop_start:envelop_stop), path.y(envelop_start:envelop_stop), 'r*');
end
%% finding contact shape

%% helper functions
function is_in_circle = getCirclePenetration(circle_x , circle_y, circle_radius, point_x , point_y)
is_in_circle = ((point_x - circle_x)^2 + (point_y - circle_y)^2) < circle_radius^2;
end

function pentration_boundaries = getPenetration(tyre_centre_x, tyre_centre_y, tyre_radius, path_x , path_y)
search_start = find(path_x > tyre_centre_x  - tyre_radius, 1 , 'first');
search_end = find(path_x > tyre_centre_x + tyre_radius , 1 , "first");
search_ind = search_start;
while search_ind < search_end
    if getCirclePenetration(tyre_centre_x , tyre_centre_y , tyre_radius , path_x(search_ind) , path_y(search_ind))
        start_ind = search_ind;
        while getCirclePenetration(tyre_centre_x , tyre_centre_y , tyre_radius , path_x(search_ind) , path_y(search_ind))
            search_ind = search_ind + 1;
        end
        stop_ind = search_ind;
        [closest_radius, closes_ind] = getClosestPoint(path_x(start_ind:stop_ind), path_y(start_ind:stop_ind),...
            tyre_centre_x, tyre_centre_y, start_ind);
        if ~exist("contact_boundaries", 'var')
            pentration_boundaries = struct("start", start_ind, "stop", stop_ind,...
                "closest_radius", closest_radius, "closest_ind", closes_ind);
        else
            pentration_boundaries(end+1) = struct("start", start_ind, "stop", stop_ind,...
                "closest_radius", closest_radius, "closest_ind", closes_ind);
        end

    end

    search_ind = search_ind+1;

end
end

function [dist , ind] = getClosestPoint(segment_x , segment_y , target_x, target_y, ind_offset)
[dist, ind] = min(sqrt((segment_x - target_x).^2 + (segment_y - target_y).^2));
ind = ind + ind_offset;
end

function [start_ind,stop_ind] = getContact(...
    path, contact_start_ind, contact_stop_ind,contact_centre_ind, angle_threshold)
start_ind = contact_centre_ind;
stop_ind = contact_centre_ind;
while(stop_ind ~= contact_stop_ind)
    stop_ind = stop_ind +1;
    if path.tangent(stop_ind) < (path.tangent(contact_centre_ind) - angle_threshold)
        %stop_ind = stop_ind - 1;
        break;
    end
end
while(start_ind ~= contact_start_ind)
    start_ind = start_ind -1;
    if path.tangent(start_ind) + pi > (path.tangent(contact_centre_ind)+pi + angle_threshold)
        %start_ind = start_ind + 1;
        break;
    end
end
end

function [circle_centre, circle_radius, start_angle, end_angle] = getDeformedSection(...
    tyre_centre_x, tyre_centre_y,tyre_radius, p_x, p_y, tangent_angle)
direction_normalt_to_contact= [-sin(tangent_angle); cos(tangent_angle)];
OP_vector = [p_x ; p_y] - [tyre_centre_x;tyre_centre_y];
circle_radius = (tyre_radius^2 - norm(OP_vector)^2)/(2*dot(OP_vector, direction_normalt_to_contact) + 2*tyre_radius);
circle_centre = [p_x; p_y] + circle_radius*direction_normalt_to_contact;
start_angle = atan2(-direction_normalt_to_contact(2) , -direction_normalt_to_contact(1));
circle_centre_local_frame = circle_centre - [tyre_centre_x;tyre_centre_y];
end_angle = atan2(circle_centre_local_frame(2) , circle_centre_local_frame(1));

end




