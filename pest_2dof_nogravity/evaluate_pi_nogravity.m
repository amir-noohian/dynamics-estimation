clc;
clear;

%% --- Load data from processed_data folder ---
data_id = 4;
filename = sprintf('./processed_data/processed_data_%d.mat', data_id);
data_file = load(filename);
data = data_file.processed_data;

theta2_v = data(1,:)';
theta4_v = data(2,:)';
thetad2_v = data(3,:)';
thetad4_v = data(4,:)';
thetadd2_v = data(5,:)';
thetadd4_v = data(6,:)';
grav2_v = data(7,:)';
grav4_v = data(8,:)';
tau2_v = data(9,:)';
tau4_v = data(10,:)';
tau2_nograv_v = tau2_v - grav2_v;
tau4_nograv_v = tau4_v - grav4_v;

g = 9.81;

N = length(theta2_v);
Y_result = zeros(2 * N, 8);
Tau_result = zeros(2 * N, 1);

for i = 1:N
    theta2 = theta2_v(i);
    theta4 = theta4_v(i);
    thetad2 = thetad2_v(i);
    thetad4 = thetad4_v(i);
    thetadd2 = thetadd2_v(i);
    thetadd4 = thetadd4_v(i);
    tau2_nograv = tau2_nograv_v(i);
    tau4_nograv = tau4_nograv_v(i);

    Y11 = 0.2e1 * cos(theta4) * thetad4 * thetad2 + cos(theta4) * thetad4 ^ 2 + 0.2e1 * sin(theta4) * thetadd2 + sin(theta4) * thetadd4;
    Y12 = 0.2e1 * cos(theta4) * thetadd2 - sin(theta4) * thetad4 ^ 2 - 0.2e1 * sin(theta4) * thetad4 * thetad2 + cos(theta4) * thetadd4;
    Y13 = thetadd4;
    Y14 = g * sin(theta4 + theta2);
    Y15 = -g * cos(theta4 + theta2);
    Y16 = thetadd2;
    Y17 = -cos(theta2) * g;
    Y18 = sin(theta2) * g;
    Y19 = tanh(20*thetad2);
    Y110 = thetad2;
    Y111 = 0;
    Y112 = 0;

    Y21 = sin(theta4) * thetadd2 - cos(theta4) * thetad2 ^ 2;
    Y22 = sin(theta4) * thetad2 ^ 2 + cos(theta4) * thetadd2;
    Y23 = thetadd2 + thetadd4;
    Y24 = g * sin(theta4 + theta2);
    Y25 = -g * cos(theta4 + theta2);
    Y26 = 0;
    Y27 = 0;
    Y28 = 0;
    Y29 = 0;
    Y210 = 0;
    Y211 = tanh(20*thetad4);
    Y212 = thetad4;

    % Y = [Y11, Y12, Y13, Y14, Y15, Y16, Y17, Y18, Y19, Y110, Y111, Y112; Y21, Y22, Y23, Y24, Y25, Y26, Y27, Y28, Y29, Y210, Y211 Y212];
    Y_nograv = [Y11, Y12, Y13, Y16, Y19, Y110, Y111, Y112; Y21, Y22, Y23, Y26, Y29, Y210, Y211 Y212];
    
    start_row = (i - 1) * 2 + 1;
    Y_result(start_row:start_row + 1, :) = Y_nograv;
    Tau_result(start_row:start_row + 1, :) = [tau2_nograv; tau4_nograv];
end

%% --- Parameter estimation ---
disp(['Condition number of Y: ', num2str(cond(Y_result))]);
[U, S, V] = svd(Y_result, 'econ');

pi_estimated_nograv = [-0.00599426383144650
-0.568073892550376
-0.359906717242286
-1.83273553525122
-2.41674809812703
-1.23437339787623
-0.269329752770356
-0.687181098898918];

Tau_pi_estimated = Y_result * pi_estimated_nograv;

Tau_error = sqrt(sum((Tau_result - Tau_pi_estimated).^2) / N);
Normal_Tau_error = Tau_error / (max(Tau_result) - min(Tau_result));
disp(["NRMSE Tau:", num2str(Normal_Tau_error)]);

%% --- Plot results ---
subplot(2,1,1);
plot(Tau_result(1:2:end));
hold on;
plot(Tau_pi_estimated(1:2:end));
legend('Recorded', 'Calculated');
title('Torque 2');

subplot(2,1,2);
plot(Tau_result(2:2:end));
hold on;
plot(Tau_pi_estimated(2:2:end));
legend('Recorded', 'Calculated');
title('Torque 4');










