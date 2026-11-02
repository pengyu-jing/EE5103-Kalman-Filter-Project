function Xh = kf_gains(A, C, Kf, Y, x0, Bu)
% KF_GAINS  Filter many runs at once with precomputed gains.
%   The gains of the Kalman filter do not depend on the data, so they are computed
%   once with kf.m and reused for every Monte Carlo run.
%   Kf : n x 1 x nk gains from kf.m,  Y : runs x nk measurements,  x0 = x^(0|-1)
%   Bu : optional constant input term (n x 1)
%   Xh : runs x nk x n, Xh(:, k+1, :) = x^(k|k)
[runs, nk] = size(Y);
n  = size(A, 1);
xm = repmat(x0(:)', runs, 1);                 % x^(k|k-1), one row per run
Xh = zeros(runs, nk, n);
for k = 1:nk
    e  = Y(:, k) - xm*C';
    xf = xm + e*Kf(:, :, k)';
    Xh(:, k, :) = reshape(xf, runs, 1, n);
    xm = xf*A';                                % x^(k+1|k) = A x^(k|k) (+ B u)
    if nargin > 5 && ~isempty(Bu)
        xm = xm + repmat(Bu(:)', runs, 1);
    end
end
end
