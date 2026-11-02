function out = part1(figdir)
% PART1  Random-walk target:  x(k+1) = x(k) + w(k),  y(k) = x(k) + v(k).
%   sigma_w = 0.1, sigma_v = 1, x(0) = 5, k = 0..N with N = 10000.
%   Produces Graphs 1-3, the bias and the variance, and the analytic checks of the report.
if nargin < 1
    figdir = 'figures';
end
N  = 10000;
sw = 0.1;  sv = 1;  x0 = 5;
A  = 1;  C = 1;  R1 = sw^2;  R2 = sv^2;
P0 = 1e5;  xhat0 = 0;                 % P(0|-1) and x^(0|-1), as in the notes
k  = 0:N;

% ---------------- simulation (one run) ----------------
rng(1);
w = sw*randn(1, N + 1);
v = sv*randn(1, N + 1);
x = zeros(1, N + 1);
x(1) = x0;
for i = 1:N
    x(i + 1) = A*x(i) + w(i);
end
y = C*x + v;

% ---------------- Kalman filter ----------------
r   = kf(A, C, R1, R2, y, xhat0, P0);
xh  = r.xf;
Pf  = reshape(r.Pf, 1, []);
Kf  = reshape(r.Kf, 1, []);
err = x - xh;

% ---------------- analytic steady state, eqs. (13)-(14) ----------------
Pbar  = (R1 + sqrt(R1^2 + 4*R1*R2))/2;      % P(k+1|k)
Kf_ss = Pbar/(Pbar + R2);                    % Kf = K (A = 1)
Pf_ss = R2*Kf_ss;                            % P(k|k)
a     = 1 - Kf_ss;                           % pole of the steady-state filter (15)

% ---------------- statistics, eq. (11) ----------------
out.bias        = mean(err);
out.variance    = mean(err.^2);
out.bias_ss     = mean(err(101:end));        % k >= 100
out.variance_ss = mean(err(101:end).^2);
out.var_meas    = mean((y - x).^2);
out.Kf_ss = Kf_ss;  out.Pf_ss = Pf_ss;  out.Pp_ss = Pbar;
out.lag   = a^2*R1/(1 - a^2);                % error budget (16)
out.noise = Kf_ss^2*R2/(1 - a^2);
out.k_conv = find(abs(Kf - Kf_ss)/Kf_ss < 0.01, 1) - 1;
out.bias_se = sqrt(Pf_ss*(1 + a)/((1 - a)*(N + 1)));

% known initial state for comparison: x^(0|-1) = 5, P(0|-1) = 0
r0  = kf(A, C, R1, R2, y, x0, 0);
Pf0 = reshape(r0.Pf, 1, []);

% ---------------- Monte Carlo check of P(k|k): 5000 runs, k = 0..60 ----------------
nmc = 5000;  kmc = 60;
rng(101);
W  = sw*randn(nmc, kmc + 1);
V  = sv*randn(nmc, kmc + 1);
X  = x0 + [zeros(nmc, 1), cumsum(W(:, 1:end - 1), 2)];
Xh = kf_gains(A, C, r.Kf(:, :, 1:kmc + 1), X + V, xhat0);
E  = X - Xh;
mse_mc  = mean(E.^2, 1);
mean_mc = mean(E, 1);

% ---------------- print ----------------
fprintf('\n=================== Part 1 ===================\n');
fprintf('analytic steady state : P(k+1|k) = %.4f, Kf = K = %.4f, P(k|k) = %.4f\n', Pbar, Kf_ss, Pf_ss);
fprintf('simulation at k = N   : Kf(N) = %.4f, P(N|N) = %.4f\n', Kf(end), Pf(end));
fprintf('Kf(0), Kf(1), Kf(2)   : %.5f, %.4f, %.4f  (within 1%% of steady state after k = %d)\n', ...
        Kf(1), Kf(2), Kf(3), out.k_conv);
fprintf('bias     k=0..N: %8.4f   k>=100: %8.4f   theory: 0 +- %.4f\n', out.bias, out.bias_ss, out.bias_se);
fprintf('variance k=0..N: %8.4f   k>=100: %8.4f   theory: %.4f\n', out.variance, out.variance_ss, Pf_ss);
fprintf('raw measurement error variance %.4f -> reduction %.1f-fold\n', out.var_meas, out.var_meas/out.variance);
fprintf('error budget (16): lag %.4f + noise %.4f = %.4f\n', out.lag, out.noise, out.lag + out.noise);
fprintf('steady-state pole %.4f, time constant %.1f samples\n', a, -1/log(a));
fprintf('Monte Carlo: max |E[e^2] - P(k|k)|/P(k|k) for k>=1: %.1f%%\n', ...
        100*max(abs(mse_mc(2:end) - Pf(2:kmc + 1))./Pf(2:kmc + 1)));

% ---------------- figures ----------------
f = figure('Name', 'Part 1 - Graph 1', 'Color', 'w', 'Position', [50 50 800 600]);
subplot(2, 1, 1);
plot(k, x, 'k-', k, xh, 'r:');
grid on;  xlim([0 N]);
xlabel('k');  ylabel('position');
legend('x(k)', 'xhat(k|k)', 'Location', 'northwest');
title('Graph 1: x(k) (solid) and xhat(k|k) (dotted)');
subplot(2, 1, 2);
kz = 0:200;
plot(kz, y(kz + 1), '.', 'Color', [0.6 0.6 0.6]);
hold on;
plot(kz, x(kz + 1), 'k-', 'LineWidth', 1);
plot(kz, xh(kz + 1), 'r:', 'LineWidth', 1.5);
hold off;  grid on;  xlim([0 200]);
xlabel('k');  ylabel('position');
legend('y(k)', 'x(k)', 'xhat(k|k)', 'Location', 'northeast');
title('Graph 1, zoom k = 0..200');
save_fig(f, figdir, 'p1_graph1');

f = figure('Name', 'Part 1 - Graphs 2 and 3', 'Color', 'w', 'Position', [50 50 900 600]);
subplot(2, 2, 1);
plot(k, Pf, 'r-', [0 N], [Pf_ss Pf_ss], 'b--');
grid on;  xlim([0 N]);  ylim([0 1.05]);
xlabel('k');  ylabel('P(k|k)');  title('Graph 2: P(k|k)');
legend('P(k|k)', sprintf('steady state %.4f', Pf_ss));
subplot(2, 2, 2);
kz = 0:60;
plot(kz, Pf(kz + 1), 'r.-', kz, Pf0(kz + 1), 'g-.', [0 60], [Pf_ss Pf_ss], 'b--');
grid on;  ylim([0 1.05]);
xlabel('k');  ylabel('P(k|k)');  title('Graph 2, zoom k = 0..60');
legend('P(0|-1) = 1e5', 'P(0|-1) = 0 (known x(0))', 'steady state');
subplot(2, 2, 3);
plot(k, Kf, 'r-', [0 N], [Kf_ss Kf_ss], 'b--');
grid on;  xlim([0 N]);  ylim([0 1.05]);
xlabel('k');  ylabel('K_f(k)');  title('Graph 3: K_f(k)');
legend('K_f(k)', sprintf('steady state %.4f', Kf_ss));
subplot(2, 2, 4);
plot(kz, Kf(kz + 1), 'r.-', [0 60], [Kf_ss Kf_ss], 'b--');
grid on;  ylim([0 1.05]);
xlabel('k');  ylabel('K_f(k)');  title('Graph 3, zoom k = 0..60');
save_fig(f, figdir, 'p1_graph23');

f = figure('Name', 'Part 1 - Monte Carlo check', 'Color', 'w', 'Position', [50 50 800 600]);
subplot(2, 1, 1);
plot(k, err, 'r-', k, 2*sqrt(Pf), 'k--', k, -2*sqrt(Pf), 'k--');
grid on;  xlim([0 N]);  ylim([-2.2 2.2]);
xlabel('k');  ylabel('x(k) - xhat(k|k)');  title('estimation error and +-2 sqrt(P(k|k))');
subplot(2, 1, 2);
kz = 0:kmc;
plot(kz, Pf(kz + 1), 'k-', kz, mse_mc, 'ro', kz, mean_mc, 'bx');
grid on;
xlabel('k');  legend('P(k|k)', 'E[(x - xhat)^2], 5000 runs', 'E[x - xhat], 5000 runs');
title('ensemble statistics versus P(k|k)');
save_fig(f, figdir, 'p1_mc');
end
