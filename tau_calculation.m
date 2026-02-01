% % %% ============================================================
% % %  Tau identification from cooling curve (SS316)
% % %  Author: ---
% % %  Purpose:
% % %   - Build tau(theta) from raw thermal decay data
% % %   - Ensure bounded, reusable tau function for fatigue post-processing
% % % ============================================================
% % close all
% % clear; clc;
% % 
% % %% ------------------ USER-DEFINED CONSTANTS ------------------
% % rho = 8000;          % Density of SS316 [kg/m^3]
% % C   = 500;           % Specific heat capacity [J/(kg.K)]
% % % SF  = 10;            % Sampling frequency [Hz]
% % window = 1;        % Moving average window [s]
% % 
% % %% ------------------ LOAD RAW THERMAL DATA -------------------
% % load 'C:\Users\mo170\OneDrive - The University of Waikato\PhD\Experiments\SS316_VAL\Tau\SS316L_260116_1343_Tau_recording\THERMAL_data.mat';
% % 
% % IRtime = IRtime(:);
% % IRtime = IRtime - IRtime(1); % Start time at zero
% % 
% % T_0D = T_0D(:);
% % T_ref_up = T_ref_up(:);
% % T_ref_bottom = T_ref_bottom(:);
% % 
% % %% ------------------ Sampling Frequency --------------------
% % SF = 1 / mode(round(diff(IRtime),9));
% % disp(['IR camera Sampling Frequency as per experimental data: ', num2str(SF),'Hz']);
% % 
% % 
% % %% --- Repair IR time vector segment-wise  ---
% % N = length(IRtime);
% % IRtime_ok = IRtime;          % start with a copy
% % 
% % if N < 2
% %     return;
% % end
% % 
% % % Find where the sequence decreases or stays the same
% % diffs = diff(IRtime);
% % bad = (diffs <= 0);
% % 
% % % Find start and end indices of every bad segment
% % % Padding with false at both ends to make edge detection easy
% % bad_padded = [false; bad; false];
% % starts = find(diff(bad_padded) == 1) - 1;   % first bad index in each run
% % ends = find(diff(bad_padded) == -1) + 1; % last bad index in each run
% % 
% % for k = 1:length(starts)
% %     seg_start = starts(k);      % index where decrease begins (1-based)
% %     seg_end = ends(k);        % last bad index in this run
% % 
% %     % The segment that needs fixing is from seg_start to seg_end (inclusive)
% %     % We interpolate using the last good point before seg_start
% %     % and the first good point after seg_end
% % 
% %     left_idx = seg_start - 1;                 % last good point before problem
% %     right_idx = seg_end + 1;                   % first good point after problem
% % 
% %     % --- Edge cases ---
% %     if left_idx < 1
% %         % Problem starts at the very beginning - forward extrapolate
% %         first_good_diff = diffs(find(~bad,1,'first'));
% %         t0 = IRtime(right_idx);
% %         IRtime_ok(1:right_idx) = t0 - first_good_diff*(right_idx:-1:0)';
% %         continue;
% %     end
% % 
% %     if right_idx > N
% %         % Problem goes to the very end - backward extrapolate
% %         last_good_diff = diffs(find(~bad,1,'last'));
% %         t0 = IRtime(left_idx);
% %         IRtime_ok(left_idx:end) = t0 + last_good_diff*(1:(N-left_idx+1))';
% %         continue;
% %     end
% % 
% %     % --- Normal case: interpolate between left_idx and right_idx ---
% %     t_left = IRtime_ok(left_idx);
% %     t_right = IRtime_ok(right_idx);
% % 
% %     % Number of points from left_idx to right_idx inclusive
% %     n_points = right_idx - left_idx + 1;
% %     new_times = linspace(t_left, t_right, n_points);
% % 
% %     IRtime_ok(left_idx:right_idx) = new_times';
% % end
% % 
% % % IR time plots
% % figure('Name', 'IR time');
% % plot(IRtime, '-k', LineWidth=0.5);
% % title('IR time',testName)
% % hold on; plot(IRtime, '.k', markersize=5);
% % hold on; plot(IRtime_ok, '.r', markersize=5);
% % legend('IRtime original','' ,'IRtime corrected','Location', 'northwest');
% % xlabel('Index'); ylabel('Time (s)')
% % grid on
% % 
% % 
% % figure('Name', 'diff (tTemp)');
% % plot(diff(IRtime), '-k', LineWidth=2);
% % title('diff (tTemp)',testName)
% % hold on; 
% % plot(diff(IRtime), '.k', markersize=5);
% % plot(diff(IRtime_ok), '.r', LineWidth=5);
% % legend('diff (IRtime)','','diff (IRtime corrected)','Location', 'northwest');
% % xlabel('Index'); ylabel('diff (tTemp)  (s)')
% % grid on
% % 
% % %% ------------------ BUILD THETA -----------------------------
% % 
% % % T_ref = 0.5 * T_ref_up + 0.5 * T_ref_bottom; % Average reference temperature
% % % T_noOffset = T_0D - T_ref; % Remove reference temperature offset
% % 
% % T_0D = movmean(T_0D, SF * window);
% % T_noOffset = T_0D - min(T_0D); 
% % 
% % % Reference temperature taken as final equilibrium temperature
% % T_inf = mean(T_noOffset(end-200:end));
% % 
% % theta = T_noOffset - T_inf;
% 
% % Enforce physical constraint
% % theta(theta < 0) = 0;
% 
% %% ------------------ SMOOTH THETA -----------------------------
% SFnew = round(window * SF);
% if mod(SFnew,2) == 0
%     SFnew = SFnew + 1;   % Ensure odd window length
% end
% 
% theta_smooth = movmean(theta, SFnew, 'endpoints','shrink');
% 
% %% ------------------ COMPUTE d(theta)/dt ---------------------
% dtheta_dt = gradient(theta_smooth, IRtime_ok); % ????????????????????????????????????????????
% 
% %% ------------------ COMPUTE TAU(theta) ----------------------
% tau_raw = nan(size(theta_smooth));
% 
% valid_idx = (theta_smooth > 0) & (dtheta_dt < 0);
% 
% tau_raw(valid_idx) = -theta_smooth(valid_idx) ./ dtheta_dt(valid_idx);
% 
% % Remove non-physical values
% tau_raw(tau_raw <= 0 | tau_raw > 1e4) = NaN;
% 
% %% ------------------ DEFINE VALID THETA RANGE ----------------
% theta_min = min(theta_smooth(valid_idx));
% theta_max = max(theta_smooth(valid_idx));
% 
% fprintf('Tau valid for theta in range: %.3f K to %.3f K\n', ...
%         theta_min, theta_max);
% 
% %% ------------------ BUILD BOUNDED TAU FUNCTION ---------------
% theta_tau_raw = theta_smooth(valid_idx);
% tau_raw_valid = tau_raw(valid_idx);
% 
% good = ~isnan(tau_raw_valid);
% theta_tau_raw = theta_tau_raw(good);
% tau_raw_valid = tau_raw_valid(good);
% 
% % ---- BINNING IN THETA SPACE ----
% dTheta = 0.05;                     % bin width [K] (adjust if needed)
% theta_edges = theta_min:dTheta:theta_max;
% theta_centres = theta_edges(1:end-1) + diff(theta_edges)/2;
% 
% tau_binned = nan(size(theta_centres));
% 
% for i = 1:length(theta_centres)
%     inBin = theta_tau_raw >= theta_edges(i) & ...
%             theta_tau_raw <  theta_edges(i+1);
%     if any(inBin)
%         tau_binned(i) = median(tau_raw_valid(inBin)); % median = robust
%     end
% end
% 
% % Remove empty bins
% good = ~isnan(tau_binned);
% theta_tau = theta_centres(good);
% tau       = tau_binned(good);
% 
% % Interpolant (bounded usage enforced later)
% tau_fun_bounded_SS316 = @(theta_q) interp1( ...
%     theta_tau, tau, ...
%     min(max(theta_q, theta_min), theta_max), ...
%     'linear', 'extrap');
% 
% %% ------------------ SAVE FOR FUTURE USE ---------------------
% save('tau_fun_bounded_SS316.mat', ...
%      'tau_fun_bounded_SS316', ...
%      'theta_min', 'theta_max', ...
%      'rho', 'C');
% 
% disp('tau_fun_bounded_SS316 saved successfully.');
% 
% %% ------------------ OPTIONAL DIAGNOSTICS --------------------
% figure;
% subplot(2,1,1)
% plot(IRtime_ok, theta_smooth, 'k');
% xlabel('Time [s]');
% ylabel('\theta [K]');
% title('Smoothed temperature decay');
% 
% subplot(2,1,2)
% plot(theta_tau, tau, '.');
% xlabel('\theta [K]');
% ylabel('\tau [s]');
% title('\tau(\theta) identified from cooling curve');
% grid on;
% 
% %%
% theta_fun = @(t,theta) -theta ./ tau_fun_bounded_SS316(theta);
% theta0 = theta_smooth(find(theta_smooth > 0,1,'first'));
% tspan  = [IRtime_ok(1), IRtime_ok(end)];
% opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
% [t_rec, theta_rec] = ode45(theta_fun, tspan, theta0, opts);
% theta_rec_interp = interp1(t_rec, theta_rec, IRtime_ok, 'linear', 'extrap');
% 
% figure;
% plot(IRtime_ok, theta_smooth, 'k', 'LineWidth',1.5); hold on;
% plot(IRtime_ok, theta_rec_interp, '--r', 'LineWidth',1.5);
% 
% xlabel('Time [s]');
% ylabel('\theta [K]');
% legend('Experimental \theta','Reconstructed \theta','Location','northeast');
% title('Verification of \tau(\theta) via temperature reconstruction');
% grid on;




%% +++++++++########################################################################################
%% ============================================================
%% ============================================================
%  Tau identification from cooling curve (SS316)
%  Author: ---
%  Purpose:
%   - Build tau(theta) from raw thermal decay data
%   - Ensure bounded, reusable tau function for fatigue post-processing
% ============================================================
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

% T_ref = 0.5 * T_ref_up + 0.5 * T_ref_bottom; % Average reference temperature
% T_noOffset = T_0D - T_ref; % Remove reference temperature offset

T_0D = movmean(T_0D, SF * window,'Endpoints','Fill');
T_inf = min(T_0D);
theta = T_0D - T_inf;

figure; plot(IRtime_ok,theta); grid on
title ('\theta')




%% Remove invalid points (remove NaN)
valid = isfinite(IRtime_ok) & isfinite(theta);
IRtime_ok     = IRtime_ok(valid);
theta = theta(valid);

%% -----------------------------
% INITIAL VALUES
% -----------------------------
theta_init = theta(1);
t_init     = IRtime_ok(1);   % usually 0

%% -----------------------------
% STEP 1: Fit theta(t)
% -----------------------------
ft = fittype('A*exp(-t/tau0) + C', ...
    'independent','t', ...
    'coefficients',{'A','tau0','C'});

A0   = theta(1) - theta(end);
tau0 = (max(IRtime_ok)-min(IRtime_ok))/5;
C0   = theta(end);

opts = fitoptions(ft);
opts.StartPoint = [A0 tau0 C0];
opts.Robust     = 'Bisquare';

[theta_fit_fun, gof] = fit(IRtime_ok, theta, ft, opts);

%% -----------------------------
% STEP 2: Reconstruct fitted theta
% -----------------------------
theta_fit = theta_fit_fun(IRtime_ok);

idx = theta_fit > 0;
Z_theta_fit = theta_fit(idx);
Z_IRtime_ok = IRtime_ok(idx);

% Convert back to temperature if needed
T_fit = Z_theta_fit + T_inf;

%% -----------------------------
% STEP 3: Calculate tau using YOUR formula
% tau(t) = -t / ln(theta/theta_init)
% -----------------------------

theta_init = Z_theta_fit(1);

ratio = Z_theta_fit / theta_init;

% Avoid log singularities
valid_tau = ratio > 0 & ratio < 1 & Z_IRtime_ok > 0;

tau_inst = NaN(size(Z_IRtime_ok));
tau_inst(valid_tau) = ...
    - Z_IRtime_ok(valid_tau) ./ log(ratio(valid_tau));

%% -----------------------------
% STEP 4: Express tau as tau(theta)
% -----------------------------
theta_tau = Z_theta_fit(valid_tau);
tau_tau   = tau_inst(valid_tau);

% Sort by theta
[theta_sorted, idx] = sort(theta_tau);
tau_sorted          = tau_tau(idx);

% Optional smoothing in theta-space
tau_theta_smooth = movmean(tau_sorted, 15);

%% -----------------------------
% STEP 5: Error metrics
% -----------------------------
residuals = theta(1:length(Z_theta_fit)) - Z_theta_fit;

RMSE  = sqrt(mean(residuals.^2));
NRMSE = RMSE / (max(theta)-min(theta));

fprintf('RMSE  = %.4f K\n', RMSE)
fprintf('NRMSE = %.4f\n', NRMSE)
fprintf('R^2   = %.4f\n', gof.rsquare)

%% -----------------------------
% PLOTS
% -----------------------------
figure
plot(IRtime_ok, theta, 'k.', 'DisplayName','Raw \theta')
hold on
plot(Z_IRtime_ok, Z_theta_fit, 'r-', 'LineWidth',1.5, ...
    'DisplayName','Fitted \theta')
xlabel('Time [s]')
ylabel('\theta [K]')
legend
grid on
title('\theta(t): raw vs fitted')

figure
plot(theta_sorted, tau_sorted, 'b.', ...
    'DisplayName','\tau(\theta)')
hold on
plot(theta_sorted, tau_theta_smooth, 'r-', 'LineWidth',1.5, ...
    'DisplayName','Smoothed \tau(\theta)')
xlabel('\theta [K]')
ylabel('\tau [s]')
legend
grid on
title('\tau as a function of \theta')
