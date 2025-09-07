clc; clear;

%% road setting
L = 3000;                             % road length
v = 15;                               % free flow speed (m/s)
dt = 1;                               % time step (s)
simulation_time = 600;                
kj = 0.15;                            % congestion density (veh/m)
w = 6;                                % shock wave speed (m/s)

kc = kj * (w / (v + w));              % cirital density
q_max = (v * w * kj) / (v + w);       % max flow (veh/s)

cell_length = v * dt;                 % cell length (m)
num_cells = floor(L / cell_length);   % cell number 
N_i = kj * cell_length;               
Q_i = q_max * dt;                     
Q_i_all = Q_i * ones(1, num_cells + 1); 

output_flow = Q_i;                    

n = zeros(simulation_time+1, num_cells);        
y = zeros(simulation_time, num_cells+1);        
R = zeros(simulation_time, num_cells+1);        

%% initial density setting
low_density_ratio = 0.45;
high_density_ratio = 0.85;

cell_0_1500   = 1 : floor(1500 / cell_length);
cell_1500_3000 = floor(1500 / cell_length)+1 : num_cells;

n(1, cell_0_1500)    = low_density_ratio  * kc * cell_length;
n(1, cell_1500_3000) = high_density_ratio * kc * cell_length;

%% main body
disturbance_cell = 100; 
disturbance_start = 150; 
disturbance_duration = 30; %T1 
recovery_period = 70; %T2

for t = 1:simulation_time
    
    % Periodic switching of input flow
    cycle_time = 200; 
    if mod(t, cycle_time) < 100
        input_ratio = 0.85;  % high density
    else
        input_ratio = 0.65;  % low density
    end
    input_flow = input_ratio * kc * cell_length;  
    
    % interference cell
    Q_i_all(:) = Q_i;  
    if t >= disturbance_start && t < disturbance_start + disturbance_duration
        Q_i_all(disturbance_cell) = 0.3 * Q_i;
    elseif t >= disturbance_start + disturbance_duration && t < disturbance_start + disturbance_duration + recovery_period
        recovery_progress = (t - (disturbance_start + disturbance_duration)) / recovery_period;
        recovered_Qi = 0.3 * Q_i + recovery_progress * ((1-0.3) * Q_i);
        Q_i_all(disturbance_cell) = recovered_Qi;
    end
    
    % Transimission
    for i = 1:num_cells+1
        if i == 1
            y(t,i) = min([input_flow, Q_i_all(i), N_i - n(t,i)]);
        elseif i == num_cells+1
            y(t,i) = min([n(t,i-1), output_flow]);
        else
            R_i = (w / v) * (N_i - n(t,i));
            R(t,i) = R_i;
            y(t,i) = min([n(t,i-1), Q_i_all(i), R_i]);
        end
    end

    % model update
    for i = 1:num_cells
        n(t+1,i) = n(t,i) + y(t,i) - y(t,i+1);
        n(t+1,i) = max(0, min(N_i, n(t+1,i)));  
    end
end

%% plot
density = n / cell_length;

figure;
imagesc(linspace(0, L, num_cells), 0:simulation_time, density);
colormap(jet); colorbar; axis xy;
xlabel('Position (m)'); ylabel('Time (s)');
title('Traffic Flow with Disturbance (85%~75%)');
clim([0 kj]);

%% congestion behaviour
congestion_threshold = kc;  
congestion_matrix = density > congestion_threshold;
congestion_length_each_t = sum(congestion_matrix, 2) * cell_length;

total_congestion_duration = sum(congestion_length_each_t > 0) * dt;
max_congestion_length = max(congestion_length_each_t);
avg_congestion_length = mean(congestion_length_each_t(congestion_length_each_t > 0));

fprintf('Total congestion duration: %.1f s\n', total_congestion_duration);
fprintf('Maximum congestion length: %.1f m\n', max_congestion_length);
fprintf('Average congestion length: %.1f m\n', avg_congestion_length);

%% plot
figure;
plot(0:simulation_time, congestion_length_each_t, 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Congested length (m)');
title('Congestion Length Over Time');
grid on;
