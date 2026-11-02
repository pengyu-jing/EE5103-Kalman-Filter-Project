function out = part3(figdir, nmc)
% PART3  Indoor positioning of a target in circular motion (25 measured samples).
%   A system of two decoupled Kalman filters (one for y, one for z) is designed:
%     M1 (baseline): constant-velocity model of Part 2, sigma_w = 0.1 and 0.7
%     M2 (final)   : harmonic-oscillator (constant-turn) model, sigma_w = 0.1,
%                    with the angular rate and centre identified from the data
%   Prints Kf(k), K(k), P(k|k), P(k+1|k) for k = 0..24, plots Figure 5.1 with the
%   estimates, and reproduces the analysis of the report (frequency response, NIS,
%   Monte Carlo).  nmc = number of Monte Carlo runs (default 5000).
if nargin < 1
    figdir = 'figures';
end
if nargin < 2
    nmc = 5000;
end
T  = 1;
R2 = 1;                         % variance of the y and z measurement noises
P0 = 1e5*eye(2);                % P(0|-1)
x0 = [0; 0];                    % x^(0|-1)
SIG_CV_SMALL = 0.1;             % M1 with the Part 2 value
SIG_CV       = 0.7;             % M1, sigma_w ~ max. centripetal acceleration r w^2
SIG_HO       = 0.1;             % M2

y = [7.1165 9.6022 8.9144 9.2717 6.3400 4.0484 0.3411 -0.6784 ...
     -5.7726 -5.4925 -9.4523 -9.7232 -9.5054 -9.7908 -7.7300 -5.9779 ...
     -4.5535 -1.5042 -0.7044 3.2406 8.3029 6.1925 9.1178 9.0904 9.0662];
z = [0.000 3.1398 6.3739 9.5877 10.1450 10.1919 9.0683 10.2254 ...
     7.5799 7.7231 5.4721 3.3990 0.9172 -1.3551 -5.2708 -9.7011 ...
     -9.4256 -9.3053 -9.3815 -9.8822 -8.1876 -8.7501 -4.5653 -1.9179 -1.0000];
n = numel(y);
k = 0:n - 1;

% ---------------- identification of the circular motion ----------------
idn = identify_circle(y, z);

% ---------------- the three designs ----------------
names  = {'M1, sigma_w = 0.1', 'M1, sigma_w = 0.7', 'M2 (oscillator)'};
models = cell(1, 3);                             % each entry: {A, G, C, R1}
[A, G, C, R1] = cv_model(T, SIG_CV_SMALL);     models{1} = {A, G, C, R1};
[A, G, C, R1] = cv_model(T, SIG_CV);           models{2} = {A, G, C, R1};
[A, G, C, R1] = osc_model(T, idn.w, SIG_HO);   models{3} = {A, G, C, R1};
centres = {[0 0], [0 0], [idn.cy idn.cz]};      % M2 filters the deviation from the centre

dof = 2*(n - 2);                                 % NIS over k = 2..24, both axes
nis_lo = 2*gammaincinv(0.025, dof/2)/dof;        % chi2 95% region of the mean NIS
nis_hi = 2*gammaincinv(0.975, dof/2)/dof;
meas_rad = hypot(y - idn.cy, z - idn.cz);

d = struct([]);
for j = 1:3
    [A, G, C, R1] = models{j}{:};
    c  = centres{j};
    ry = kf(A, C, R1, R2, y - c(1), x0, P0);
    rz = kf(A, C, R1, R2, z - c(2), x0, P0);
    d(j).ry = ry;
    d(j).ey = ry.xf(1, :) + c(1);                   % x^_y(k|k)
    d(j).ez = rz.xf(1, :) + c(2);                   % x^_z(k|k)
    d(j).ss = kf_steady(A, C, R1, R2);
    [d(j).H, d(j).h2] = kf_freq(A, C, d(j).ss, idn.w);
    nis = [ry.e(3:end).^2 ./ reshape(ry.S(1, 1, 3:end), 1, []), ...
           rz.e(3:end).^2 ./ reshape(rz.S(1, 1, 3:end), 1, [])];
    d(j).nis = mean(nis);
    rad = hypot(d(j).ey - idn.cy, d(j).ez - idn.cz);
    d(j).rad_mean = mean(rad(6:end));               % k >= 5
    d(j).rad_rms  = sqrt(mean((rad(6:end) - idn.r).^2));
    d(j).rmse_th  = sqrt(2*R2*d(j).h2 + idn.r^2*abs(d(j).H - 1)^2);   % eq. (25)
end

% ---------------- Monte Carlo on a synthetic circle with the identified geometry ----------------
ty = idn.cy + idn.r*cos(idn.th0 + idn.w*k);
tz = idn.cz + idn.r*sin(idn.th0 + idn.w*k);
rng(3);
Yn = repmat(ty, nmc, 1) + sqrt(R2)*randn(nmc, n);
Zn = repmat(tz, nmc, 1) + sqrt(R2)*randn(nmc, n);
se = zeros(6, n);                                 % raw, M1(0.1), M1(0.7), M2, M2 with 1.1w, M2 oracle
se(1, :) = mean((Yn - ty).^2 + (Zn - tz).^2, 1);
for j = 1:2
    [A, G, C, R1] = models{j}{:};
    g  = kf(A, C, R1, R2, zeros(1, n), x0, P0);   % gains do not depend on the data
    Ey = kf_gains(A, C, g.Kf, Yn, x0);
    Ez = kf_gains(A, C, g.Kf, Zn, x0);
    se(j + 1, :) = mean((Ey(:, :, 1) - ty).^2 + (Ez(:, :, 1) - tz).^2, 1);
end
[Ao, Go, Co, R1o] = osc_model(T, idn.w, SIG_HO);
go = kf(Ao, Co, R1o, R2, zeros(1, n), x0, P0);
Ey = kf_gains(Ao, Co, go.Kf, Yn - idn.cy, x0);
Ez = kf_gains(Ao, Co, go.Kf, Zn - idn.cz, x0);
se(6, :) = mean((Ey(:, :, 1) + idn.cy - ty).^2 + (Ez(:, :, 1) + idn.cz - tz).^2, 1);
id_runs = zeros(nmc, 3);
for i = 1:nmc
    idr = identify_circle(Yn(i, :), Zn(i, :));    % identification repeated on every run
    id_runs(i, :) = [idr.cy, idr.cz, idr.w];
    for m = 1:2
        wm = idr.w*(1 + 0.1*(m - 1));
        [A, G, C, R1] = osc_model(T, wm, SIG_HO);
        ry = kf(A, C, R1, R2, Yn(i, :) - idr.cy, x0, P0);
        rz = kf(A, C, R1, R2, Zn(i, :) - idr.cz, x0, P0);
        se(3 + m, :) = se(3 + m, :) + ((ry.xf(1, :) + idr.cy - ty).^2 + (rz.xf(1, :) + idr.cz - tz).^2)/nmc;
    end
end
rmse   = sqrt(se);
rmse_m = sqrt(mean(se(:, 6:end), 2));             % k = 5..24

% ---------------- print ----------------
fprintf('\n=================== Part 3 ===================\n');
fprintf('identified circle: centre (%.2f, %.2f), radius %.2f\n', idn.cy, idn.cz, idn.r);
fprintf('angular rate %.4f +- %.4f rad/sample (%.2f deg/sample, period %.2f samples)\n', ...
        idn.w, idn.se_w, idn.w*180/pi, 2*pi/idn.w);
fprintf('centripetal acceleration r w^2 = %.2f\n', idn.r*idn.w^2);
fprintf('A_w = [%.4f %.4f; %.4f %.4f]\n', models{3}{1}');
fprintf('\n%-30s %10s %10s %10s %8s\n', '', 'M1 (0.1)', 'M1 (0.7)', 'M2', 'raw');
fprintf('%-30s %10.3f %10.3f %10.3f %8.3f\n', 'P11 (Riccati)', d(1).ss.Pf(1, 1), d(2).ss.Pf(1, 1), d(3).ss.Pf(1, 1), R2);
fprintf('%-30s %10.3f %10.3f %10.3f %8.3f\n', 'R2 ||H||_2^2', R2*d(1).h2, R2*d(2).h2, R2*d(3).h2, R2);
fprintf('%-30s %10.3f %10.3f %10.3f %8.3f\n', '|H(e^jw)|', abs(d(1).H), abs(d(2).H), abs(d(3).H), 1);
fprintf('%-30s %10.1f %10.1f %10.1f %8.1f\n', 'angle H(e^jw) (deg)', angle(d(1).H)*180/pi, ...
        angle(d(2).H)*180/pi, angle(d(3).H)*180/pi, 0);
fprintf('%-30s %10.2f %10.2f %10.2f %8.2f\n', 'radius r|H|', idn.r*abs(d(1).H), idn.r*abs(d(2).H), ...
        idn.r*abs(d(3).H), idn.r);
fprintf('%-30s %10.2f %10.2f %10.2f %8.2f\n', 'mean radius of estimates', d(1).rad_mean, d(2).rad_mean, ...
        d(3).rad_mean, mean(meas_rad(6:end)));
fprintf('%-30s %10.2f %10.2f %10.2f %8.2f\n', 'RMS radial deviation', d(1).rad_rms, d(2).rad_rms, ...
        d(3).rad_rms, sqrt(mean((meas_rad(6:end) - idn.r).^2)));
fprintf('%-30s %10.2f %10.2f %10.2f %8s\n', sprintf('mean NIS [%.2f, %.2f]', nis_lo, nis_hi), ...
        d(1).nis, d(2).nis, d(3).nis, '-');
fprintf('%-30s %10.2f %10.2f %10.2f %8.2f\n', 'MC RMSE, k = 5..24', rmse_m(2), rmse_m(3), rmse_m(4), rmse_m(1));
fprintf('%-30s %10.2f %10.2f %10.2f %8.2f\n', 'MC RMSE, k = 24', rmse(2, end), rmse(3, end), rmse(4, end), rmse(1, end));
fprintf('%-30s %10.2f %10.2f %10.2f %8.2f\n', 'RMSE theory, eq. (25)', d(1).rmse_th, d(2).rmse_th, ...
        d(3).rmse_th, sqrt(2*R2));
fprintf('M2 with 1.1 w: RMSE %.2f;  M2 with the true parameters (k = 24): %.2f\n', rmse_m(5), rmse(6, end));
fprintf('MC standard errors of the identification: cy %.2f, cz %.2f, w %.4f\n', std(id_runs));

print_gain_table(1, d(3), 'Final design M2 (oscillator, sigma_w = 0.1): identical for the y- and z-filters');
print_gain_table(1, d(2), 'Baseline M1 (constant velocity, sigma_w = 0.7): identical for the y- and z-filters');
fprintf('\nEstimates of the final design M2\n%4s %10s %10s %10s %10s\n', 'k', 'y(k)', 'z(k)', 'xhat_y', 'xhat_z');
fprintf('%4d %10.4f %10.4f %10.4f %10.4f\n', [k; y; z; d(3).ey; d(3).ez]);

% the same tables as text files
resdir = fullfile(fileparts(figdir), 'results');
if ~exist(resdir, 'dir')
    mkdir(resdir);
end
fid = fopen(fullfile(resdir, 'part3_tables.txt'), 'w');
print_gain_table(fid, d(3), 'Final design M2 (oscillator, sigma_w = 0.1): identical for the y- and z-filters');
print_gain_table(fid, d(2), 'Baseline M1 (constant velocity, sigma_w = 0.7): identical for the y- and z-filters');
fprintf(fid, '\nEstimates of the final design M2\n%4s %10s %10s %10s %10s\n', 'k', 'y(k)', 'z(k)', 'xhat_y', 'xhat_z');
fprintf(fid, '%4d %10.4f %10.4f %10.4f %10.4f\n', [k; y; z; d(3).ey; d(3).ez]);
fclose(fid);

% ---------------- self-check against the report (deterministic part) ----------------
ref = [0.2615583624779451, 0.3290021478026292, 0.2773347967554649, 10.241015505874469, ...
       9.819548007710724, 0.0354303079316135, 0.2786785891672717, 0.2784825666150622, ...
       0.04508871527152529, 0.18682345955181773];
val = [idn.w, idn.cy, idn.cz, idn.r, d(3).ey(end), d(3).ez(end), d(3).ry.Pf(1, 1, end), ...
       d(3).ss.Kf(1), d(3).ss.Kf(2), d(3).h2];
if max(abs(val - ref)./abs(ref)) < 1e-6
    fprintf('\nSelf-check against the report: PASS\n');
else
    fprintf('\nSelf-check against the report: FAIL (max rel. deviation %.2e)\n', max(abs(val - ref)./abs(ref)));
end

out.idn = idn;  out.designs = d;  out.rmse = rmse;  out.rmse_mean = rmse_m;

% ---------------- figures ----------------
f = figure('Name', 'Part 3 - Figure 5.1 with estimates', 'Color', 'w', 'Position', [50 50 600 600]);
plot(y, z, 'o', 'Color', [0 0.447 0.741]);
hold on;
plot(d(3).ey, d(3).ez, 'rx-', 'LineWidth', 1, 'MarkerSize', 8);
hold off;
axis([-15 15 -15 15]);  axis square;  grid on;
xlabel('y-axis');  ylabel('z-axis');
legend('measurement (y, z)', 'estimate (xhat_y(k|k), xhat_z(k|k))', 'Location', 'northwest');
title('Figure 5.1 with the estimates of the final design (M2)');
save_fig(f, figdir, 'p3_fig51');

f = figure('Name', 'Part 3 - comparison of the designs', 'Color', 'w', 'Position', [50 50 1200 420]);
t = linspace(0, 2*pi, 361);
for j = 1:3
    subplot(1, 3, j);
    plot(idn.cy + idn.r*cos(t), idn.cz + idn.r*sin(t), '--', 'Color', [0.6 0.6 0.6]);
    hold on;
    plot(y, z, 'o', 'Color', [0 0.447 0.741]);
    plot(d(j).ey, d(j).ez, 'rx-');
    hold off;
    axis([-15 15 -15 15]);  axis square;  grid on;
    xlabel('y-axis');  ylabel('z-axis');  title(names{j});
end
save_fig(f, figdir, 'p3_compare');

f = figure('Name', 'Part 3 - gains and covariances', 'Color', 'w', 'Position', [50 50 1000 600]);
labels = {'K_{f1}(k)', 'K_{f2}(k)', 'P_{11}(k|k)', 'P_{22}(k|k)'};
for i = 1:4
    subplot(2, 2, i);
    vals = zeros(2, n);
    sel  = [3 2];                                         % M2 and M1 (sigma_w = 0.7)
    for m = 1:2
        r = d(sel(m)).ry;
        switch i
            case 1, vals(m, :) = reshape(r.Kf(1, 1, :), 1, []);
            case 2, vals(m, :) = reshape(r.Kf(2, 1, :), 1, []);
            case 3, vals(m, :) = reshape(r.Pf(1, 1, :), 1, []);
            case 4, vals(m, :) = reshape(r.Pf(2, 2, :), 1, []);
        end
    end
    plot(k, vals(1, :), 'r.-', k, vals(2, :), 'bs--', 'MarkerSize', 4);
    grid on;  xlim([0 24]);  xlabel('k');  title(labels{i});
    if i == 4
        set(gca, 'YScale', 'log');
    end
    if i == 1
        legend('M2 (oscillator)', 'M1 (CV), sigma_w = 0.7');
    end
end
save_fig(f, figdir, 'p3_gains');

f = figure('Name', 'Part 3 - Monte Carlo RMSE', 'Color', 'w', 'Position', [50 50 800 420]);
plot(k, rmse(1, :), 'k-', k, rmse(2, :), 'g--', k, rmse(3, :), 'bs-', k, rmse(4, :), 'ro-', ...
     k, rmse(5, :), 'r:', 'LineWidth', 1);
grid on;  xlim([0 24]);  xlabel('k');  ylabel('position RMSE');
legend('raw measurement', 'M1, sigma_w = 0.1', 'M1, sigma_w = 0.7', 'M2, identified w', 'M2, 1.1 w', ...
       'Location', 'northwest');
title(sprintf('position RMSE over %d synthetic runs', nmc));
save_fig(f, figdir, 'p3_mc');
end


function print_gain_table(fid, dj, caption)
% Kf(k), K(k), P(k|k), P(k+1|k) for k = 0..24 and the steady state.
r = dj.ry;
fprintf(fid, '\n%s\n', caption);
fprintf(fid, '%4s %9s %9s %9s %9s | %9s %9s %11s | %11s %11s %11s\n', 'k', 'Kf1', 'Kf2', 'K1', 'K2', ...
        'P11(k|k)', 'P12(k|k)', 'P22(k|k)', 'P11(k+1|k)', 'P12(k+1|k)', 'P22(k+1|k)');
for i = 1:size(r.Kf, 3)
    fprintf(fid, '%4d %9.4f %9.4f %9.4f %9.4f | %9.4f %9.4f %11.4g | %11.4g %11.4g %11.4g\n', i - 1, ...
            r.Kf(1, 1, i), r.Kf(2, 1, i), r.K(1, 1, i), r.K(2, 1, i), ...
            r.Pf(1, 1, i), r.Pf(1, 2, i), r.Pf(2, 2, i), ...
            r.Pp(1, 1, i + 1), r.Pp(1, 2, i + 1), r.Pp(2, 2, i + 1));
end
s = dj.ss;
fprintf(fid, '%4s %9.4f %9.4f %9.4f %9.4f | %9.4f %9.4f %11.4g | %11.4g %11.4g %11.4g\n', 'inf', ...
        s.Kf(1), s.Kf(2), s.K(1), s.K(2), s.Pf(1, 1), s.Pf(1, 2), s.Pf(2, 2), s.Pp(1, 1), s.Pp(1, 2), s.Pp(2, 2));
end
