function out = part2(figdir)
% PART2  Constant-velocity target, eqs. (5.3)-(5.4):
%   x(k+1) = [1 T; 0 1] x(k) + [T^2/2; T] w(k),  y(k) = [1 0] x(k) + v(k)
%   sigma_w = 0.1, sigma_v = 1, T = 1, x(0) = [0; 30], k = 0..N with N = 10000.
%   Produces Graphs 4-9, the biases and variances, and the analytic checks of the report.
if nargin < 1
    figdir = 'figures';
end
N  = 10000;
T  = 1;  sw = 0.1;  sv = 1;
[A, G, C, R1] = cv_model(T, sw);
R2 = sv^2;
x0 = [0; 30];
P0 = 1e5*eye(2);  xhat0 = [0; 0];
k  = 0:N;

% ---------------- simulation (one run) ----------------
rng(2);
w = sw*randn(1, N + 1);
v = sv*randn(1, N + 1);
x = zeros(2, N + 1);
x(:, 1) = x0;
for i = 1:N
    x(:, i + 1) = A*x(:, i) + G*w(i);
end
y = C*x + v;

% ---------------- Kalman filter ----------------
r   = kf(A, C, R1, R2, y, xhat0, P0);
xh  = r.xf;
err = x - xh;
P11 = reshape(r.Pf(1, 1, :), 1, []);
P22 = reshape(r.Pf(2, 2, :), 1, []);
Kf1 = reshape(r.Kf(1, 1, :), 1, []);
Kf2 = reshape(r.Kf(2, 1, :), 1, []);

% ---------------- analytic steady state (alpha-beta filter), eqs. (19)-(22) ----------------
lam   = sw*T^2/sv;                                   % tracking index
u     = (4 + lam - sqrt(lam^2 + 8*lam))/4;
alpha = 1 - u^2;
beta  = lam*u;
Pf_an = sv^2*[alpha, beta/T; beta/T, beta*(alpha - beta/2)/((1 - alpha)*T^2)];
S     = sv^2/(1 - alpha);
Pp_an = [alpha*S, beta*S/T; beta*S/T, Pf_an(2, 2) + sw^2*T^2];
Kf_an = [alpha; beta/T];
K_an  = A*Kf_an;
ss    = kf_steady(A, C, R1, R2);                     % numerical check of the closed form

% ---------------- statistics, eq. (11) ----------------
out.bias        = mean(err, 2);
out.variance    = mean(err.^2, 2);
out.bias_ss     = mean(err(:, 101:end), 2);          % k >= 100
out.variance_ss = mean(err(:, 101:end).^2, 2);
out.first_term  = err(2, 1)^2/(N + 1);               % contribution of k = 0 to the velocity variance
out.alpha = alpha;  out.beta = beta;  out.Pf_an = Pf_an;  out.Pp_an = Pp_an;
out.k_conv = find(all(abs([Kf1; Kf2] - Kf_an)./Kf_an < 0.01, 1), 1) - 1;
poles = eig(A - K_an*C);
vdiff = y(2:end) - y(1:end - 1);                     % naive velocity estimate
var_vdiff = mean((x(2, 1:end - 1) - vdiff).^2);

% ---------------- Monte Carlo check of P(k|k): 5000 runs, k = 0..40 ----------------
nmc = 5000;  kmc = 40;
rng(102);
W  = sw*randn(nmc, kmc + 1);
V  = sv*randn(nmc, kmc + 1);
Xs = zeros(nmc, kmc + 1, 2);
Xs(:, 1, 1) = x0(1);
Xs(:, 1, 2) = x0(2);
for i = 1:kmc
    Xs(:, i + 1, 1) = Xs(:, i, 1) + T*Xs(:, i, 2) + G(1)*W(:, i);
    Xs(:, i + 1, 2) = Xs(:, i, 2) + G(2)*W(:, i);
end
Xh  = kf_gains(A, C, r.Kf(:, :, 1:kmc + 1), Xs(:, :, 1) + V, xhat0);
mse = squeeze(mean((Xs - Xh).^2, 1));                % (kmc+1) x 2

% ---------------- print ----------------
fprintf('\n=================== Part 2 ===================\n');
fprintf('tracking index lambda = %.2f -> u = %.2f, alpha = %.4f, beta = %.4f\n', lam, u, alpha, beta);
fprintf('P(k|k)   analytic = [%.4f %.4f; %.4f %.4f], Riccati iteration = [%.4f %.4f; %.4f %.4f]\n', ...
        Pf_an(1, 1), Pf_an(1, 2), Pf_an(2, 1), Pf_an(2, 2), ss.Pf(1, 1), ss.Pf(1, 2), ss.Pf(2, 1), ss.Pf(2, 2));
fprintf('P(k+1|k) analytic = [%.4f %.4f; %.4f %.4f]\n', Pp_an(1, 1), Pp_an(1, 2), Pp_an(2, 1), Pp_an(2, 2));
fprintf('Kf = [%.4f; %.4f], K = A Kf = [%.4f; %.4f]\n', Kf_an(1), Kf_an(2), K_an(1), K_an(2));
fprintf('simulation at k = N: P11 = %.4f, P12 = %.4f, P22 = %.4f, Kf1 = %.4f, Kf2 = %.4f\n', ...
        r.Pf(1, 1, end), r.Pf(1, 2, end), r.Pf(2, 2, end), Kf1(end), Kf2(end));
fprintf('Kf(0) = [%.4f; %.4f], Kf(1) = [%.4f; %.4f], P22(1|1) = %.4f; within 1%% after k = %d\n', ...
        Kf1(1), Kf2(1), Kf1(2), Kf2(2), P22(2), out.k_conv);
fprintf('poles of A - K C: %.4f +- j%.4f (modulus %.4f)\n', real(poles(1)), abs(imag(poles(1))), abs(poles(1)));
fprintf('bias      x1: %8.4f (k>=100: %8.4f)   x2: %8.4f (k>=100: %8.4f)\n', ...
        out.bias(1), out.bias_ss(1), out.bias(2), out.bias_ss(2));
fprintf('variance  x1: %8.4f (k>=100: %8.4f)   x2: %8.4f (k>=100: %8.4f)\n', ...
        out.variance(1), out.variance_ss(1), out.variance(2), out.variance_ss(2));
fprintf('k = 0 term of the velocity variance: %.4f; naive differencing variance %.3f\n', ...
        out.first_term, var_vdiff);
fprintf('Monte Carlo: max relative deviation from P11, P22 for k>=1: %.1f%%, %.1f%%\n', ...
        100*max(abs(mse(2:end, 1)' - P11(2:kmc + 1))./P11(2:kmc + 1)), ...
        100*max(abs(mse(2:end, 2)' - P22(2:kmc + 1))./P22(2:kmc + 1)));

% ---------------- figures ----------------
f = figure('Name', 'Part 2 - Graphs 4 and 5', 'Color', 'w', 'Position', [50 50 1000 650]);
subplot(2, 2, 1);
plot(k, x(1, :), 'k-', k, xh(1, :), 'r:');
grid on;  xlim([0 N]);
xlabel('k');  ylabel('position');  title('Graph 4: x_1(k) (solid), xhat_1(k|k) (dotted)');
legend('x_1(k)', 'xhat_1(k|k)', 'Location', 'northwest');
subplot(2, 2, 3);
plot(k, err(1, :), 'r-', k, 2*sqrt(P11), 'k--', k, -2*sqrt(P11), 'k--');
grid on;  xlim([0 N]);  ylim([-3 3]);
xlabel('k');  ylabel('x_1 - xhat_1');  title('Graph 4 detail: position error and +-2 sqrt(P_{11})');
subplot(2, 2, 2);
plot(k, x(2, :), 'k-', k, xh(2, :), 'r:');
grid on;  xlim([0 N]);
xlabel('k');  ylabel('velocity');  title('Graph 5: x_2(k) (solid), xhat_2(k|k) (dotted)');
legend('x_2(k)', 'xhat_2(k|k)', 'Location', 'northwest');
subplot(2, 2, 4);
kz = 0:60;
plot(kz, x(2, kz + 1), 'k-', kz, xh(2, kz + 1), 'r:.', 'LineWidth', 1.2);
grid on;  xlim([0 60]);
xlabel('k');  ylabel('velocity');  title('Graph 5, zoom k = 0..60');
save_fig(f, figdir, 'p2_graph45');

f = figure('Name', 'Part 2 - Graphs 6 to 9', 'Color', 'w', 'Position', [50 50 1000 800]);
items = {P11, Pf_an(1, 1), 'P_{11}(k|k)', 'Graph 6'; P22, Pf_an(2, 2), 'P_{22}(k|k)', 'Graph 7'; ...
         Kf1, Kf_an(1), 'K_{f1}(k)', 'Graph 8';    Kf2, Kf_an(2), 'K_{f2}(k)', 'Graph 9'};
kz = 0:30;
for i = 1:4
    val = items{i, 1};  ssv = items{i, 2};
    subplot(4, 2, 2*i - 1);
    plot(k, val, 'r-', [0 N], [ssv ssv], 'b--');
    grid on;  xlim([0 N]);
    ylabel(items{i, 3});  title(sprintf('%s: %s (steady state %.4f)', items{i, 4}, items{i, 3}, ssv));
    subplot(4, 2, 2*i);
    plot(kz, val(kz + 1), 'r.-', [0 30], [ssv ssv], 'b--');
    grid on;  xlim([0 30]);
    title(sprintf('%s, zoom k = 0..30', items{i, 4}));
    if i == 2
        set(gca, 'YScale', 'log');
        subplot(4, 2, 2*i - 1);
        set(gca, 'YScale', 'log');
    end
end
subplot(4, 2, 7);  xlabel('k');
subplot(4, 2, 8);  xlabel('k');
save_fig(f, figdir, 'p2_graph6789');

f = figure('Name', 'Part 2 - Monte Carlo check', 'Color', 'w', 'Position', [50 50 800 400]);
kz = 0:kmc;
semilogy(kz, P11(kz + 1), 'k-', kz, P22(kz + 1), 'k--', kz, mse(:, 1), 'ro', kz, mse(:, 2), 'bs');
grid on;  xlabel('k');  ylabel('variance');
legend('P_{11}(k|k)', 'P_{22}(k|k)', 'E[(x_1-xhat_1)^2]', 'E[(x_2-xhat_2)^2]');
title('diagonal of P(k|k) versus mean-square errors over 5000 runs');
save_fig(f, figdir, 'p2_mc');
end
