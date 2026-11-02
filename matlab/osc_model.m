function [A, G, C, R1] = osc_model(T, w, sigma_w)
% OSC_MODEL  Harmonic-oscillator (constant turn-rate) model of one axis of a circular motion.
%   State [p - c; dp/dt], exact discretisation of  d2p/dt2 = -w^2 (p - c),  eq. (24) of the report.
th = w*T;
A  = [cos(th), sin(th)/w; -w*sin(th), cos(th)];
G  = [T^2/2; T];
C  = [1 0];
R1 = sigma_w^2*(G*G');
end
