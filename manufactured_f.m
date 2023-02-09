function [fx, fy] = manufactured_f(x,y,t,nu)

    u = t * (1-cos(2*pi*x)).*y.*(2-3.*y);
    v = -t * 2* pi * sin(2*pi*x).*(y.^2).*(1-y);

    fx = (1 - (y.*(2-3*y)).*cos(2*pi*x)) + u.* (2*pi*t*sin(2*pi*x)).*y.*(2-3*y) + v.*t*(1-cos(2*pi*x)).*(2-6*y) + 2*nu*(x-0.5) - nu * (4*pi^2*t*cos(2*pi*x).*(y.*(2-3*y)) - 6*t*(1-cos(2*pi*x)));

    fy = -2*pi*sin(2*pi*x).*(y.^2).*(1-y) + u.*(-t*4*pi^2*cos(2*pi*x).*(y.^2).*(1-y)) + v.*(-t*2*pi*sin(2*pi*x).*(2*y-3*y.^2)) + 2*nu*(y-0.5) - nu*(t*8*pi^3*sin(2*pi*x).*(y.^2).*(1-y) - t*2*pi*sin(2*pi*x).*(2-6*y));

    fx = reshape(fx, numel(fx),1);
    fy = reshape(fy, numel(fy),1);

end