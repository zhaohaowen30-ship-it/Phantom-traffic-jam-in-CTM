clc; clear; close all;

% [Kp, Ki, Kd, lower_limit, upper_limit]
lb = [0, 0, 0, 0.3, 0.7];   % Lower Bound
ub = [1, 0.1, 0.2, 0.7, 1.0];   % uppon bound
x0 = [0.5, 0.05, 0.1, 0.5, 1.0];  % initial value

options = optimoptions('ga','MaxGenerations',20,'PopulationSize',40,'Display','iter');
best_params = ga(@(x)obj_fun(x),5,[],[],[],[],lb,ub,[],options);

fprintf('\nresult:\n Kp=%.3f, Ki=%.3f, Kd=%.3f, lower=%.3f, upper=%.3f\n', ...
    best_params(1), best_params(2), best_params(3), best_params(4), best_params(5));

%% finial result
[total_duration, max_length, avg_length, density, congestion_length_each_t, L, num_cells, simulation_time, kj] = ...
    run_simulation(best_params(1), best_params(2), best_params(3), best_params(4), best_params(5));

fprintf('\nfinial result:\n');
fprintf('  Times: %.1f s\n', total_duration);
fprintf('  Max length: %.1f m\n', max_length);
fprintf('  Avg length: %.1f m\n', avg_length);

%% plot
figure;
imagesc(linspace(0, L, num_cells), 0:simulation_time, density);
colormap(jet); colorbar; axis xy;
xlabel('position (m)'); ylabel('time (s)');
title('traffic density (final optimized)');
clim([0 kj]);

[~, ~, ~, ~, congestion_length_no_ctrl, ~, ~, ~, ~] = ...
    run_simulation(0, 0, 0, 1, 1);   % without control

figure;
plot(0:simulation_time, congestion_length_no_ctrl, 'r-', 'LineWidth', 1.5); hold on;
plot(0:simulation_time, congestion_length_each_t, 'b-', 'LineWidth', 1.5);
xlabel('time (s)'); ylabel('congestion length (m)');
title('Congestion length over time: Control vs No Control');
legend('No Control','With Control');
grid on;



%% cost function
function cost = obj_fun(params)

    Kp = params(1); Ki = params(2); Kd = params(3);
    lower_limit = params(4); upper_limit = params(5);

    [total_duration, max_length, avg_length] = run_simulation(Kp, Ki, Kd, lower_limit, upper_limit);

    w1 = 1/3; w2 = 1/3; w3 = 1/3;
    cost = w1*total_duration + w2*max_length + w3*avg_length;
end


%% 
function [total_congestion_duration, max_congestion_length, avg_congestion_length, density, congestion_length_each_t, L, num_cells, simulation_time, kj] ...
    = run_simulation(Kp, Ki, Kd, lower_limit, upper_limit)

    %% basic setting
    L = 3000;                             % road length (m)
    v = 15;                               % free flow speed (m/s)
    dt = 1;                               % time step (s)
    simulation_time = 600;                % total time  (s)
    kj = 0.15;                            % congestion density (veh/m)
    w = 6;                                % shock wave speed (m/s)

    kc = kj * (w / (v + w));              % cirtal density
    q_max = (v * w * kj) / (v + w);       % max flow (veh/s)

    cell_length = v * dt;                 % cell length (m)
    num_cells = floor(L / cell_length);   % cell number
    N_i = kj * cell_length;              
    Q_i = q_max * dt;                     
    Q_i_all = Q_i * ones(1, num_cells + 1); 

    output_flow = Q_i;                    % output flow


    n = zeros(simulation_time+1, num_cells);

    %% inital density
    initial_density_ratio = 0.3;
    n(1,:) = initial_density_ratio * kc * cell_length;

    %% main part
    disturbance_cell = 100; 
    disturbance_start = 150; 
    disturbance_duration = 30; 
    recovery_period = 70; 

    setpoint = 1;   
    lock_cell = 0;
    integral_error = 0;
    prev_error = 0;

    for t = 1:simulation_time
     
    cycle_time = 200;  
    if mod(t, cycle_time) < 100
        input_ratio = 0.85;  % high density
    else
        input_ratio = 0.75;  % low density
    end
    input_flow = input_ratio * kc * cell_length; %input flow 

    % interference cell setting
    Q_i_all(:) = Q_i;  

        if t >= disturbance_start && t < disturbance_start + disturbance_duration
            Q_i_all(disturbance_cell) = 0.3 * Q_i;
        elseif t >= disturbance_start + disturbance_duration && ...
               t < disturbance_start + disturbance_duration + recovery_period
            recovery_progress = (t - (disturbance_start + disturbance_duration)) / recovery_period;
            recovered_Qi = 0.3 * Q_i + recovery_progress * ((1-0.3) * Q_i); 
            Q_i_all(disturbance_cell) = recovered_Qi;
        end

        % Call the PID controller
        [Q_i_all, lock_cell, integral_error, prev_error] = pid_controller(n(t,:), Q_i_all, Q_i, ...
            Kp, Ki, Kd, setpoint, lock_cell, integral_error, prev_error, lower_limit, upper_limit);

        % transimission between cells
        y = zeros(1,num_cells+1);
        for i = 1:num_cells+1
            if i == 1
                y(i) = min([input_flow, Q_i_all(i), N_i - n(t,i)]);
            elseif i == num_cells+1
                y(i) = min([n(t,i-1), output_flow]);
            else
                R_i = (w / v) * (N_i - n(t,i));
                y(i) = min([n(t,i-1), Q_i_all(i), R_i]);
            end
        end

        % update
        for i = 1:num_cells
            n(t+1,i) = n(t,i) + y(i) - y(i+1);
            n(t+1,i) = max(0, min(N_i, n(t+1,i)));  
        end
    end

    %% congestion 
    density = n / cell_length;
    congestion_threshold = kc;  
    congestion_matrix = density > congestion_threshold;
    congestion_length_each_t = sum(congestion_matrix, 2) * cell_length;

    total_congestion_duration = sum(congestion_length_each_t > 0) * dt;
    max_congestion_length = max(congestion_length_each_t);
    avg_congestion_length = mean(congestion_length_each_t(congestion_length_each_t > 0));
end


%% PID function
function [Q_i_all_updated, lock_cell_updated, integral_error_updated, prev_error_updated] = ...
    pid_controller(n_t, Q_i_all, Q_i, Kp, Ki, Kd, setpoint, lock_cell, integral_error, prev_error, lower_limit, upper_limit)

    num_cells = length(n_t);
    dt = 1;   

    Q_i_all_updated = Q_i_all;
    lock_cell_updated = lock_cell;
    integral_error_updated = integral_error;
    prev_error_updated = prev_error;

    if lock_cell == 0
        for i = 1:num_cells
            if n_t(i) > setpoint
                lock_cell_updated = i;
                integral_error_updated = 0;
                prev_error_updated = 0;
                break
            end
        end
    else

        current_value = n_t(lock_cell);
        error = (current_value - setpoint)/setpoint;

        integral_error_updated = integral_error + error * dt;  
        derivative_error = (error - prev_error) / dt;          
        prev_error_updated = error;

        adjustment = Kp * error + Ki * integral_error_updated + Kd * derivative_error;

        if adjustment > 0
            new_coefficient = max(lower_limit, 1 - adjustment);
        else
            new_coefficient = min(upper_limit, 1 + adjustment);
        end

        start_idx = max(lock_cell - 60, 1);
        end_idx = lock_cell - 1;
        if end_idx >= start_idx
            Q_i_all_updated(start_idx:end_idx) = new_coefficient * Q_i;
        end

        if current_value < 0.8 * setpoint
            lock_cell_updated = 0;
        end
    end
end
