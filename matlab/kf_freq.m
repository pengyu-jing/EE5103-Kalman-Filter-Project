function [H, h2] = kf_freq(A, C, ss, w)
% KF_FREQ  Steady-state filter from y to x^_1(k|k), eq. (10) of the report.
%   H  = H(exp(j w)),  h2 = squared H2 norm of H (sum of the squared impulse response).
n = size(A, 1);
z = exp(1i*w);
Hm = (eye(n) - ss.Kf*C) * ((z*eye(n) - A + ss.K*C) \ ss.K) + ss.Kf;
H  = Hm(1, 1);
% impulse response: h(0) = Kf1, h(i) = Co F^(i-1) K with F = A - K C, Co = row 1 of (I - Kf C)
F  = A - ss.K*C;
Co = eye(n) - ss.Kf*C;
Co = Co(1, :);
h2 = ss.Kf(1)^2;
v  = ss.K;
for i = 1:5000
    h2 = h2 + (Co*v)^2;
    v  = F*v;
end
end
