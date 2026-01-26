%% ============================================================
%  Tau identification from cooling curve (SS316)
%  Author: ---
%  Purpose:
%   - Build tau(theta) from raw thermal decay data
%   - Ensure bounded, reusable tau function for fatigue post-processing
% ============================================================

clear; clc;

%% ------------------ USER-DEFINED CONSTANTS ------------------
rho = 8000;          % Density of SS316 [kg/m^3]
C   = 500;           % Specific heat capacity [J/(kg.K)]
% SF  = 10;            % Sampling frequency [Hz]
window = 1;        % Moving average window [s]

%% ------------------ LOAD RAW THERMAL DATA -------------------
load THERMAL_data.mat;

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

T_ref = 0.5 * T_ref_up + 0.5 * T_ref_bottom; % Average reference temperature
T_noOffset = T_0D - T_ref; % Remove reference temperature offset

% Reference temperature taken as final equilibrium temperature
T_inf = mean(T_noOffset(end-200:end));

theta = T_noOffset - T_inf;

% Enforce physical constraint
theta(theta < 0) = 0;

%% ------------------ SMOOTH THETA -----------------------------
SFnew = round(window * SF);
if mod(SFnew,2) == 0
    SFnew = SFnew + 1;   % Ensure odd window length
end

theta_smooth = movmean(theta, SFnew, 'endpoints','shrink');

%% ------------------ COMPUTE d(theta)/dt ---------------------
dtheta_dt = gradient(theta_smooth, IRtime_ok); % ????????????????????????????????????????????

%% ------------------ COMPUTE TAU(theta) ----------------------
tau_raw = nan(size(theta_smooth));

valid_idx = (theta_smooth > 0) & (dtheta_dt < 0);

tau_raw(valid_idx) = -theta_smooth(valid_idx) ./ dtheta_dt(valid_idx);

% Remove non-physical values
tau_raw(tau_raw <= 0 | tau_raw > 1e4) = NaN;

%% ------------------ DEFINE VALID THETA RANGE ----------------
theta_min = min(theta_smooth(valid_idx));
theta_max = max(theta_smooth(valid_idx));

fprintf('Tau valid for theta in range: %.3f K to %.3f K\n', ...
        theta_min, theta_max);

%% ------------------ BUILD BOUNDED TAU FUNCTION ---------------
theta_tau = theta_smooth(valid_idx);
tau   = tau_raw(valid_idx);

% Remove NaNs
good = ~isnan(tau);
theta_tau = theta_tau(good);
tau   = tau(good);

% Sort (required for interpolation)
[theta_tau, sortIdx] = sort(theta_tau);
tau = tau(sortIdx);

% Interpolant (bounded usage enforced later)
tau_fun_bounded_SS316 = @(theta_q) interp1( ...
    theta_tau, tau, ...
    min(max(theta_q, theta_min), theta_max), ...
    'linear', 'extrap');

%% ------------------ SAVE FOR FUTURE USE ---------------------
save('tau_fun_bounded_SS316.mat', ...
     'tau_fun_bounded_SS316', ...
     'theta_min', 'theta_max', ...
     'rho', 'C');

disp('tau_fun_bounded_SS316 saved successfully.');

%% ------------------ OPTIONAL DIAGNOSTICS --------------------
figure;
subplot(2,1,1)
plot(IRtime_ok, theta_smooth, 'k');
xlabel('Time [s]');
ylabel('\theta [K]');
title('Smoothed temperature decay');

subplot(2,1,2)
plot(theta_tau, tau, '.');
xlabel('\theta [K]');
ylabel('\tau [s]');
title('\tau(\theta) identified from cooling curve');
grid on;
