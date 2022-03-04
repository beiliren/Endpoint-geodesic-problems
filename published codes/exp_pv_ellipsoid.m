% exp_pv_ellipsoid.m
function r1 = exp_pv_ellipsoid(x,v,opts)

A = opts.A;
eps0 = opts.eps0;
C_S = opts.C_S;

d = size(x,1);
if sqrt(v'*A*v) <= eps0
    r1 = x;
else
    [~,sol] = ode45(@geodesiceq,[0,1],[x;v]);
    r1 = sol(end,1:d)';
end

    function drdt = geodesiceq(t,r)
        % geodesic equation
        drdt = [r(d+1:2*d);-C_S(r(1:d),r(d+1:2*d),r(d+1:2*d))];
    end

end