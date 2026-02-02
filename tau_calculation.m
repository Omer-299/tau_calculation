
close all
clear; clc;

%% ------------------ USER-DEFINED CONSTANTS ------------------
rho = 8000;          % Density of SS316 [kg/m^3]
C   = 500;           % Specific heat capacity [J/(kg.K)]
window = 3;        % Moving average window [s]

%% ------------------ LOAD RAW THERMAL DATA -------------------
load 'C:\Users\mo170\OneDrive - The University of Waikato\PhD\Experiments\SS316_VAL\Tau\SS316L_260116_1343_Tau_recording\THERMAL_data.mat';

IRtime = IRtime(:);
IRtime = IRtime - IRtime(1); % Start time at zero

T_0D = T_0D(:);
T_ref_up = T_ref_up(:);
T_ref_bottom = T_ref_bottom(:);

%% ------------------ Sampling Frequency --------------------
SF = 1 / mode(round(diff(IRtime),9));
disp(['IR camera Sampling Frequency as per experimental data: ', num2str(SF),'Hz']);


%% --- Repair IR time vector segment-wise  ---
N = length(IRtime);
IRtime_ok = IRtime;          % start with a copy

if N < 2
    return;
end

% Find where the sequence decreases or stays the same
diffs = diff(IRtime);
bad = (diffs <= 0);

% Find start and end indices of every bad segment
% Padding with false at both ends to make edge detection easy
bad_padded = [false; bad; false];
starts = find(diff(bad_padded) == 1) - 1;   % first bad index in each run
ends = find(diff(bad_padded) == -1) + 1; % last bad index in each run

for k = 1:length(starts)
    seg_start = starts(k);      % index where decrease begins (1-based)
    seg_end = ends(k);        % last bad index in this run

    % The segment that needs fixing is from seg_start to seg_end (inclusive)
    % We interpolate using the last good point before seg_start
    % and the first good point after seg_end

    left_idx = seg_start - 1;                 % last good point before problem
    right_idx = seg_end + 1;                   % first good point after problem

    % --- Edge cases ---
    if left_idx < 1
        % Problem starts at the very beginning - forward extrapolate
        first_good_diff = diffs(find(~bad,1,'first'));
        t0 = IRtime(right_idx);
        IRtime_ok(1:right_idx) = t0 - first_good_diff*(right_idx:-1:0)';
        continue;
    end

    if right_idx > N
        % Problem goes to the very end - backward extrapolate
        last_good_diff = diffs(find(~bad,1,'last'));
        t0 = IRtime(left_idx);
        IRtime_ok(left_idx:end) = t0 + last_good_diff*(1:(N-left_idx+1))';
        continue;
    end

    % --- Normal case: interpolate between left_idx and right_idx ---
    t_left = IRtime_ok(left_idx);
    t_right = IRtime_ok(right_idx);

    % Number of points from left_idx to right_idx inclusive
    n_points = right_idx - left_idx + 1;
    new_times = linspace(t_left, t_right, n_points);

    IRtime_ok(left_idx:right_idx) = new_times';
end

% IR time plots
figure('Name', 'IR time');
plot(IRtime, '-k', LineWidth=0.5);
title('IR time',testName)
hold on; plot(IRtime, '.k', markersize=5);
hold on; plot(IRtime_ok, '.r', markersize=5);
legend('IRtime original','' ,'IRtime corrected','Location', 'northwest');
xlabel('Index'); ylabel('Time (s)')
grid on


figure('Name', 'diff (tTemp)');
plot(diff(IRtime), '-k', LineWidth=2);
title('diff (tTemp)',testName)
hold on; 
plot(diff(IRtime), '.k', markersize=5);
plot(diff(IRtime_ok), '.r', LineWidth=5);
legend('diff (IRtime)','','diff (IRtime corrected)','Location', 'northwest');
xlabel('Index'); ylabel('diff (tTemp)  (s)')
grid on

%% ------------------ BUILD THETA -----------------------------

T_inf = min(T_0D);
diff_Tt_Tinf = T_0D - T_inf;

figure; plot(IRtime_ok,diff_Tt_Tinf); grid on
title ('\theta')

%% -----------------------------
% INITIAL VALUES
% -----------------------------
diff_Tinit_Tinf = T_0D(1) - T_inf;
t_init     = IRtime_ok(1);   % usually 0

%% -----------------------------
% STEP 1: Fit theta(t)
% -----------------------------
ft = fittype('A*exp(-t/tau0) + C', ...
    'independent','t', ...
    'coefficients',{'A','tau0','C'});

A0   = diff_Tt_Tinf(1) - diff_Tt_Tinf(end);
tau0 = (max(IRtime_ok)-min(IRtime_ok))/5;
C0   = diff_Tt_Tinf(end);

opts = fitoptions(ft);
opts.StartPoint = [A0 tau0 C0];
opts.Robust     = 'Bisquare';

[theta_fit_fun, gof] = fit(IRtime_ok, diff_Tt_Tinf, ft, opts);

%% -----------------------------
% STEP 2: Reconstruct fitted theta
% -----------------------------
theta_fit = theta_fit_fun(IRtime_ok);

idx = theta_fit > 0;
theta_fit = theta_fit(idx);
IRtime_fit = IRtime_ok(idx);

% Convert back to temperature if needed
T_fit = theta_fit + T_inf;

%% -----------------------------
% STEP 3: Calculate tau using YOUR formula
% tau(t) = -t / ln(theta/theta_init)
% -----------------------------

diff_Tinit_Tinf = theta_fit(1);

ratio = theta_fit / diff_Tinit_Tinf;

% Avoid log singularities
valid_tau = ratio > 0 & ratio < 1 & IRtime_fit > 0;

tau_inst = NaN(size(IRtime_fit));
tau_inst(valid_tau) = ...
    - IRtime_fit(valid_tau) ./ log(ratio(valid_tau));

%% -----------------------------
% STEP 4: Express tau as tau(theta)
% -----------------------------
theta_tau = theta_fit(valid_tau);
tau   = tau_inst(valid_tau);

%% -----------------------------
% STEP 5: Error metrics
% -----------------------------
residuals = diff_Tt_Tinf(1:length(theta_fit)) - theta_fit;

RMSE  = sqrt(mean(residuals.^2));
NRMSE = RMSE / (max(diff_Tt_Tinf)-min(diff_Tt_Tinf));

fprintf('RMSE  = %.4f K\n', RMSE)
fprintf('NRMSE = %.4f\n', NRMSE)
fprintf('R^2   = %.4f\n', gof.rsquare)

%% -----------------------------
% PLOTS
% -----------------------------
figure
plot(IRtime_ok, diff_Tt_Tinf, 'k.', 'DisplayName','Raw \theta')
hold on
plot(IRtime_fit, theta_fit, 'r-', 'LineWidth',1.5, ...
    'DisplayName','Fitted \theta')
xlabel('Time [s]')
ylabel('\theta [K]')
legend
grid on
title('\theta(t): raw vs fitted')

%%
figure
plot(theta_tau, tau, 'b.', ...
    'DisplayName','\tau(\theta)')
xlabel('\theta [K]')
ylabel('\tau [s]')
legend
grid on
title('\tau as a function of \theta')


%% ================ Saving Tau(theta)
tau_table = NaN(length(tau),2);
tau_table(:,1) = theta_tau;
tau_table(:,2) = tau;

save('tau_table.mat','tau_table')
