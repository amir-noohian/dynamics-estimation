function pi_nogravity(data_id, pi_id, mode)
% Usage:
%   estimate_or_evaluate_nograv(4, 4, 'calc')  % Use data 4, save pi as 4
%   estimate_or_evaluate_nograv(4, 10, 'calc') % Use data 4, save pi as 10
%   estimate_or_evaluate_nograv(4, 10, 'eval') % Use data 4, load pi 10

clc;
fprintf("Loading processed data for data_id = %d\n", data_id);

%% Load processed data
filename = sprintf('./processed_data/processed_data_%d.mat', data_id);
if ~isfile(filename)
    error('Processed data file not found: %s', filename);
end

data_file = load(filename);
data = data_file.processed_data;

theta2_v    = data(1,:)';
theta4_v    = data(2,:)';
thetad2_v   = data(3,:)';
thetad4_v   = data(4,:)';
thetadd2_v  = data(5,:)';
thetadd4_v  = data(6,:)';
grav2_v     = data(7,:)';
grav4_v     = data(8,:)';
tau2_v      = data(9,:)';
tau4_v      = data(10,:)';

% Remove gravity from torque
tau2_nograv_v = tau2_v - grav2_v;
tau4_nograv_v = tau4_v - grav4_v;

N = length(theta2_v);
Y_result = zeros(2 * N, 8);
Tau_result = zeros(2 * N, 1);

% Build regressor matrix
for i = 1:N
    theta2   = theta2_v(i);
    theta4   = theta4_v(i);
    thetad2  = thetad2_v(i);
    thetad4  = thetad4_v(i);
    thetadd2 = thetadd2_v(i);
    thetadd4 = thetadd4_v(i);
    tau2_ng  = tau2_nograv_v(i);
    tau4_ng  = tau4_nograv_v(i);

    Y11 = 2 * cos(theta4) * thetad4 * thetad2 + cos(theta4) * thetad4^2 + 2 * sin(theta4) * thetadd2 + sin(theta4) * thetadd4;
    Y12 = 2 * cos(theta4) * thetadd2 - sin(theta4) * thetad4^2 - 2 * sin(theta4) * thetad4 * thetad2 + cos(theta4) * thetadd4;
    Y13 = thetadd4;
    Y14 = thetadd2;
    Y15 = tanh(20 * thetad2);
    Y16 = thetad2;
    Y17 = 0;
    Y18 = 0;

    Y21 = sin(theta4) * thetadd2 - cos(theta4) * thetad2^2;
    Y22 = sin(theta4) * thetad2^2 + cos(theta4) * thetadd2;
    Y23 = thetadd2 + thetadd4;
    Y24 = 0;
    Y25 = 0;
    Y26 = 0;
    Y27 = tanh(20 * thetad4);
    Y28 = thetad4;

    Y = [Y11, Y12, Y13, Y14, Y15, Y16, Y17, Y18;
         Y21, Y22, Y23, Y24, Y25, Y26, Y27, Y28];

    idx = (i - 1) * 2 + 1;
    Y_result(idx:idx+1, :) = Y;
    Tau_result(idx:idx+1)  = [tau2_ng; tau4_ng];
end

%% Calculate or load pi_estimated_nograv
pi_path = sprintf('./processed_data/pi_estimated_nograv_%d.mat', pi_id);

switch lower(mode)
    case 'calc'
        fprintf("Calculating new pi_estimated_nograv for pi_id = %d...\n", pi_id);
        [U, S, V] = svd(Y_result, 'econ');
        pi_estimated_nograv = V * ((U' * Tau_result) ./ diag(S));

        save(pi_path, 'pi_estimated_nograv');
        fprintf("Saved pi_estimated_nograv to %s\n", pi_path);

    case 'eval'
        fprintf("Loading pi_estimated_nograv from %s\n", pi_path);
        if ~isfile(pi_path)
            error("File does not exist. Run in 'calc' mode first.");
        end
        pi_struct = load(pi_path);
        pi_estimated_nograv = pi_struct.pi_estimated_nograv;

    otherwise
        error("Mode must be 'calc' or 'eval'.");
end

%% Evaluate torque prediction
Tau_pi_estimated = Y_result * pi_estimated_nograv;
Tau_error = sqrt(sum((Tau_result - Tau_pi_estimated).^2) / N);
Normal_Tau_error = Tau_error / (max(Tau_result) - min(Tau_result));

fprintf("NRMSE Torque Error: %.6f\n", Normal_Tau_error);
fprintf("pi_estimated_nograv:\n");
disp(pi_estimated_nograv);

%% Plot results
subplot(2,1,1);
plot(Tau_result(1:2:end)); hold on;
plot(Tau_pi_estimated(1:2:end), '--');
legend('Measured', 'Estimated');
title(sprintf('Torque Joint 2 without Gravity (Data ID: %d, PI ID: %d)', data_id, pi_id));

subplot(2,1,2);
plot(Tau_result(2:2:end)); hold on;
plot(Tau_pi_estimated(2:2:end), '--');
legend('Measured', 'Estimated');
title(sprintf('Torque Joint 4 without Gravity(Data ID: %d, PI ID: %d)', data_id, pi_id));

end
