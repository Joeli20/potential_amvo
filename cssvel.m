function [u,w] = cssvel(x,z,l)
sigma = 1;

r1 = sqrt(x^2 + z^2);

r2 = sqrt((x-l)^2 + z^2);

omega1 = atan(x/z);
omega2 = atan((x-l)/z);

u = sigma/(4*pi) * log(r1^2/r2^2);

w = sigma/(2*pi) * (omega2 - omega1);

end