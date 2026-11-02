function r = kf(A, C, R1, R2, y, x0, P0, Bu)
% KF  Kalman filter in the notation of the EE5103 notes, eqs. (1.52)-(1.57).
%
%   r = kf(A, C, R1, R2, y, x0, P0)      runs the filter on y(0),...,y(N)
%   r = kf(A, C, R1, R2, y, x0, P0, Bu)  with a known input term B*u(k)
%
%   y  : p x (N+1) measurements, column k+1 holds y(k)
%   x0 : x^(0|-1),  P0 : P(0|-1)
%   Bu : n x (N+1) known input term B*u(k) (optional)
%
%   Column / page k+1 of every output corresponds to time k:
%   r.xf = x^(k|k)   r.xp = x^(k|k-1)  (N+2 columns, last one is x^(N+1|N))
%   r.Pf = P(k|k)    r.Pp = P(k|k-1)   (N+2 pages)
%   r.Kf, r.K        gains (n x p x (N+1))
%   r.e  = y(k) - C x^(k|k-1) (innovation),  r.S = C P(k|k-1) C' + R2
n  = size(A, 1);
p  = size(C, 1);
nk = size(y, 2);
if nargin < 8 || isempty(Bu)
    Bu = zeros(n, nk);
end
r.xf = zeros(n, nk);       r.xp = zeros(n, nk + 1);
r.Pf = zeros(n, n, nk);    r.Pp = zeros(n, n, nk + 1);
r.Kf = zeros(n, p, nk);    r.K  = zeros(n, p, nk);
r.e  = zeros(p, nk);       r.S  = zeros(p, p, nk);
r.xp(:, 1)    = x0(:);
r.Pp(:, :, 1) = P0;
for k = 1:nk
    Pm = r.Pp(:, :, k);                          % P(k|k-1)
    S  = C*Pm*C' + R2;
    Kf = Pm*C'/S;                                % (1.52)
    K  = A*Pm*C'/S;                              % (1.53)
    e  = y(:, k) - C*r.xp(:, k);                 % innovation
    r.xf(:, k)     = r.xp(:, k) + Kf*e;          % (1.54)  x^(k|k)
    r.xp(:, k + 1) = A*r.xp(:, k) + Bu(:, k) + K*e;   % (1.55)  x^(k+1|k)
    P  = Pm - Pm*C'/S*C*Pm;                      % (1.56)  P(k|k)
    Pn = A*Pm*A' - K*S*K' + R1;                  % (1.57)  P(k+1|k)
    r.Pf(:, :, k)     = (P + P')/2;
    r.Pp(:, :, k + 1) = (Pn + Pn')/2;
    r.Kf(:, :, k) = Kf;
    r.K(:, :, k)  = K;
    r.e(:, k)     = e;
    r.S(:, :, k)  = S;
end
end
