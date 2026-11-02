function idn = identify_circle(y, z)
% IDENTIFY_CIRCLE  Least-squares identification of a uniform circular motion.
%   (i)  Kasa circle fit: min sum (y^2 + z^2 + D y + E z + F)^2, linear in (D, E, F)
%   (ii) angular rate: straight-line fit of the unwrapped bearing, theta(k) = th0 + w k
y = y(:);
z = z(:);
k = (0:numel(y) - 1)';
sol = [y, z, ones(size(y))] \ (-(y.^2 + z.^2));
idn.cy = -sol(1)/2;
idn.cz = -sol(2)/2;
idn.r  = sqrt(idn.cy^2 + idn.cz^2 - sol(3));
th = unwrap(atan2(z - idn.cz, y - idn.cy));
p  = polyfit(k, th, 1);
idn.w   = p(1);
idn.th0 = p(2);
res = th - polyval(p, k);
n = numel(k);
idn.se_w = sqrt(sum(res.^2)/(n - 2)) / sqrt(sum((k - mean(k)).^2));   % standard error of w
end
