function pro3_muchchunoor_subhashchandra(path_planning_strategy_type)
    % -------------------------
    % CS5331: Aerial Computing
    % Project #3 
    % Summer II, 2023
    % -------------------------

    r = 100; % Coverage range of circle radius

    % Network size: width (w) x height (h)
    w_begin = 0;
    w_end = 1000;
    h_begin = 0;
    h_end = 1000;

    if path_planning_strategy_type<0 || path_planning_strategy_type>3
        % Terminate the function
        fprintf('Invalid input, choose the proper input as per the below instructions:\n0 for Random Planning \n1 for Distance-based Planning\n2 for Density-based Planning\n3 for Clustered Graph\n');
        return;
    end
    % Number of cells (subareas): 5-by-5, by default
    n_cell = 5;
    tot_cell = n_cell * n_cell;
    size_cell = w_end / n_cell;

    % Number of target points: 
    n = 100;

    % Set of rectangular subareas: 5-by-5
    x_ = linspace(w_begin, w_end - size_cell, n_cell);
    ux = [];
    for i = 1:n_cell
        ux = [ux, x_]; 
    end 
    ux = ux';

    y_ = ones(1, n_cell);  
    uy = [];
    for i = 1:n_cell
        uy = [uy, y_ .* (size_cell * (i - 1))];
    end 
    uy = uy';

    % Number of weights: w, n-by-1, uniform
    w_begin = 0;
    w_end = 10;
    w = w_begin + (w_end - w_begin) .* rand(n, 1);

    % -----------------------
    % Scale-free distribution 
    % -----------------------

    % Clustering exponent, alpha
    alpha = 1.4;

    % Population, pop  
    % -- Initialize to zero
    pop = ones(tot_cell, 1) - 1;

    % Probability, prob
    % -- Initialize to zero
    prob = ones(tot_cell, 1) - 1;

    % A set of rectangular subareas, 25-by-5
    subarea_loc = [ux, uy];

    % The first target point is randomly assigned to one of cells
    pos_subarea = randi(tot_cell);
    pos_x = randi(size_cell) + ux(pos_subarea);
    pos_y = randi(size_cell) + uy(pos_subarea);
    pop(pos_subarea) = pop(pos_subarea) + 1;

    % The first target point - randomly assigned
    loc(1, 1) = pos_x;
    loc(1, 2) = pos_y;

    % Generate all scale-free target points (x, y)
    for i = 2:n
        % Calculate probabilities
        sigma_pop = 0;
        for j = 1: tot_cell
            sigma_pop = sigma_pop + power(pop(j) + 1, alpha);
        end
        for j = 1: tot_cell
            prob(j) = power(pop(j) + 1, alpha) / sigma_pop;
        end

        % Choose one of subareas based on the probability
        rand_prob = rand(1, 1); % Generate between 0 to 1
        cumu_prob = 0; 
        for j = 1: tot_cell
            cumu_prob = cumu_prob + prob(j);
            if (cumu_prob >= rand_prob)
                pos_subarea = j;
                break
            end
        end

        % Generate a position within the chosen subarea
        pos_x = randi(size_cell) + ux(pos_subarea);
        pos_y = randi(size_cell) + uy(pos_subarea);
        % Increment the population of subarea
        pop(pos_subarea) = pop(pos_subarea) + 1;

        % Add a new target point's (x, y) into a row
        loc = [loc; [pos_x, pos_y]];
    end

    % Do the Neighbor Density-based Clustering on the scale-free points
    clustered_points = neighbor_density_clustering(loc, r);

    % Get the number of POIs (n) and their locations (loc)
    n_poi = size(loc, 1);

    % filter the scan points from clustered points
    scan_points = cell2mat(clustered_points);

    % Plot target points (POIs), scan points, and coverage circles
    plot(loc(:, 1), loc(:, 2), "x", 'MarkerSize', 10, 'LineWidth', 1.5,'Color','g')
    hold on
    plot(scan_points(:, 1), scan_points(:, 2), "^", 'MarkerSize', 8, 'LineWidth', 2.0)
    for i = 1:size(scan_points, 1)
        x_circ = scan_points(i, 1);
        y_circ = scan_points(i, 2);
        theta = linspace(0, 2 * pi, 100);
        x_circle = r * cos(theta) + x_circ;
        y_circle = r * sin(theta) + y_circ;
        plot(x_circle, y_circle, 'r--', 'LineWidth', 1);
    end

    % Mark a star symbol for the base point (0,0)
    plot(0, 0, '*', 'MarkerSize', 12, 'Color', 'red','LineWidth', 2.5,'DisplayName', 'base(0,0)')

    set(gca, 'FontSize', 16);
    set(gca, 'xTick', [-200:200:1400]);
    set(gca, 'yTick', [-200:200:1400]);
    xlabel('X', 'FontSize', 16);
    ylabel('Y', 'FontSize', 16); 
    axis([-200 1400 -200 1400]);

    legend('POIs', 'scan points', 'coverage range','FontSize', 10, 'Location', 'northeast')
    title('After clustering and before path planning strategy', 'FontSize', 12);
    axis equal;

    if path_planning_strategy_type == 3
        % Exit the function
        fprintf('Graph plotted after doing the Neighbor Density Clustering algorithm for 100 POIs on 1000x1000 m^2 network.\n');
        return;
    end

    % Initialize the drone start position at (0, 0)
    current_x = 0;
    current_y = 0;

    % Initialize the sequence of visits (path) for the drone with the starting position
    path = [current_x, current_y];
    
    % different path planning strategies
    while ~isempty(scan_points)
        if path_planning_strategy_type == 0   % Random Planning (RAND)
            next_idx = randi(size(scan_points, 1));
            path_planning_title='Random Planning (RAND)';
        elseif path_planning_strategy_type == 1  % Nearest Neighbor First (NNF)
            distances = sqrt((scan_points(:, 1) - current_x).^2 + (scan_points(:, 2) - current_y).^2);
            [~, next_idx] = min(distances);
            path_planning_title='Nearest Neighbor First (NNF)';
        elseif path_planning_strategy_type == 2  % Density-based First (DF)
            covered_pois = count_poi_in_coverage_radius(scan_points(:, 1), scan_points(:, 2), loc(:, 1), loc(:, 2), r);
            poi_counts = sum(covered_pois, 1);
            [~, next_idx] = max(poi_counts);
            path_planning_title='Density-based First (DF)';
        end

        % Update the current position and add the selected scan point to the path
        current_x = scan_points(next_idx, 1);
        current_y = scan_points(next_idx, 2);
        path = [path; [current_x, current_y]];
        scan_points(next_idx, :) = [];
    end

    % Plot drone path on the same graph
    plot(path(:, 1), path(:, 2), 'Color','b', 'LineWidth', 2,'DisplayName', 'drone path')
    hold off
    title(strcat(path_planning_title,' path planning strategy'), 'FontSize', 12);
    % drone flying distance from the base(0,0) to the last scan point
    total_distance = sum(sqrt(diff(path(:, 1)).^2 + diff(path(:, 2)).^2));

    fprintf('The followed path planning strategy is: %s\n', path_planning_title);

    fprintf('Total flying distance of the Drone from the network base(0,0): %.2f (m)\n', total_distance);
end

%{ 
This function will take the POIs and coverage area radius as inputs and-
-perform the clustering and then returns clustered points.
%}
function clustered_points = neighbor_density_clustering(target_points, coverage_area_radius)
    neighbor_density = zeros(size(target_points, 1), 1);
    for i = 1:size(target_points, 1)
        for j = 1:size(target_points, 1)
            if i ~= j
                d = sqrt((target_points(i,1) - target_points(j,1))^2 + (target_points(i,2) - target_points(j,2))^2);
                if d <= coverage_area_radius
                    neighbor_density(i) = neighbor_density(i) + 1;
                end
            end
        end
    end

    [~, sorting_idx] = sort(neighbor_density);
    sorted_target_points = target_points(sorting_idx, :);
    clustered_points = {};
    distance_cluster_points = distance_based_clustering(sorted_target_points, coverage_area_radius);
    clustered_points = [clustered_points, distance_cluster_points];
end

%{
This method will do the distance clustering which will be used in the-
-neighbor density based clustering algorithm.
%}
function cluster_points = distance_based_clustering(target_points, coverage_area_radius)
    cluster_points = [];
    c_init = target_points(randi(size(target_points, 1)), :);
    cluster_points = [cluster_points; c_init];
    for i = 1:size(target_points, 1)
        xi = target_points(i, 1);
        yi = target_points(i, 2);
        covered = false;
        for j = 1:size(cluster_points, 1)
            xj = cluster_points(j, 1);
            yj = cluster_points(j, 2);
            d = sqrt((xi - xj)^2 + (yi - yj)^2);
            if d <= coverage_area_radius
                covered = true;
                break;
            end
        end

        if ~covered
            cluster_points = [cluster_points; [xi, yi]];
        end
    end
end

%{
This function is to count the number of POIs in the coverage area radius-
- and it is used in the Density-based First path planning strategy.
%}
function num_poi = count_poi_in_coverage_radius(x_circle, y_circle, x_points, y_points, radius)
    distances = sqrt((x_points - x_circle').^2 + (y_points - y_circle').^2);
    num_poi = sum(distances <= radius);
end
