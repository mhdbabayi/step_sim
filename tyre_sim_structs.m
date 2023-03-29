constants_s = struct(...
    "beta", 1, ...
    "separation_threshold",deg2rad(10), ...
    "tyre_num_nodes", 360,...
    "tyre_radius", 1, ...
    "road_total_length", 5, ...
    "step_position", 3, ...
    "step_width", 0.3, ...
    "step_height", 0.1, ...
    "step_phase_width", pi, ... 
    "road_distance_step", 0.01 ...
    );

tyre_s = struct(...
    "free_radius", constants_s.tyre_radius,...
    "theta", linspace(0 , 2*pi * (1 - 1/constants_s.tyre_num_nodes), constants_s.tyre_num_nodes)',...
    "r_free", zeros(constants_s.tyre_num_nodes , 1) + constants_s.tyre_radius, ...
    "r_deformed", zeros(constants_s.tyre_num_nodes, 1) + constants_s.tyre_radius,...
    "x_centre", constants_s.step_position - 0.3 ,...
    "y_centre", constants_s.tyre_radius+0.01, ...
    "penetration_inds", [] ,... indices in vector theta  corresponding to penetrations
    "contact_centre_inds", [] , ... indices in vector theta corresponding to the closest point in each contact to the tyre centre
    "contact_boundary_inds", [] ... enveloping
    );
road_s = struct(...
    "x", [0:constants_s.road_distance_step:constants_s.road_total_length]', ...
    "y", [0:constants_s.road_distance_step:constants_s.road_total_length]' * 0,...
    "gradient", [0:constants_s.road_distance_step:constants_s.road_total_length]' * 0,...
    "ddy", [],...
    "curvature",[0:constants_s.road_distance_step:constants_s.road_total_length]' * 0,  ...
    "penetration_inds", [], ...
    "contact_centre_inds", [] , ... indices in x and y corresponding to the closest point in each contact to the tyre centre
    "contact_boundary_inds", [] ... enveloping
    ) ;
plot_handles_s = struct(...
    "frame" , [], ...
    "static", []...
    );
%% for auto completion
tyre = tyre_s;
road = road_s;
constants = constants_s;
plot_handles = plot_handles_s;