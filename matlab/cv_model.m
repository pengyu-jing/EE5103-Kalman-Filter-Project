function [A, G, C, R1] = cv_model(T, sigma_w)
% CV_MODEL  Constant-velocity (white-noise acceleration) model, eqs. (5.3)-(5.4).
%   State [position; velocity],  x(k+1) = A x(k) + G w(k),  y(k) = C x(k) + v(k).
A  = [1 T; 0 1];
G  = [T^2/2; T];
C  = [1 0];
R1 = sigma_w^2*(G*G');
end
