function ss = kf_steady(A, C, R1, R2)
% KF_STEADY  Steady-state solution of the Riccati equation (1.57).
%   Iterates (1.57) until convergence, so no toolbox is needed (the result is the
%   same as idare/dare of the Control System Toolbox).
%   ss.Pp = P(k+1|k), ss.Pf = P(k|k), ss.Kf, ss.K at steady state.
n  = size(A, 1);
Pm = eye(n);
for it = 1:100000
    S  = C*Pm*C' + R2;
    K  = A*Pm*C'/S;
    Pn = A*Pm*A' - K*S*K' + R1;
    Pn = (Pn + Pn')/2;
    done = max(abs(Pn(:) - Pm(:))) < 1e-15*max(1, max(abs(Pn(:))));
    Pm = Pn;
    if done
        break
    end
end
S     = C*Pm*C' + R2;
ss.Pp = Pm;
ss.Kf = Pm*C'/S;
ss.K  = A*ss.Kf;
ss.Pf = Pm - ss.Kf*C*Pm;
end
