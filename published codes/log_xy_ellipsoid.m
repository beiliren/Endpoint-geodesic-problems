% log_xy_ellipsoid.m
function solout = log_xy_ellipsoid(x,y,T,options)


d = size(x,1);
if sqrt((x-y)'*options.A*(x-y)) < options.eps0
    solout.x = 0;
    solout.y = [x;options.TPj(x,y,T)];
else
    dt = linspace(0,T,21);
    solinit = bvpinit(dt,@ini_guess);
    opts = bvpset('FJacobian',@jac);%,'RelTol',0.1,'AbsTol',0.1,'NMax',1000
    solout = bvp4c(@bvpfcn,@bcfcn,solinit,opts); % bvp5c
end

    function g = ini_guess(tt)
        g = [x;options.TPj(x,y,T)*tt];
    end

    function drdt = bvpfcn(t,r)
        % geodesic equation
        drdt = [r(d+1:2*d);-options.C_S(r(1:d),r(d+1:2*d),r(d+1:2*d))];
    end

    function res = bcfcn(r0,r1)
        % boundary conditions
        res = [r0(1:d)-x;r1(1:d)-y];
    end

    function dfdt = jac(t,r)
        dfdt = [zeros(d,d),eye(d);
            -options.DC_S(r(1:d),[1;0;0],r(d+1:2*d)),-options.DC_S(r(1:d),[0;1;0],r(d+1:2*d)),-options.DC_S(r(1:d),[0;0;1],r(d+1:2*d)),-2*options.C_S(r(1:d),r(d+1:2*d),[1;0;0]),-2*options.C_S(r(1:d),r(d+1:2*d),[0;1;0]),-2*options.C_S(r(1:d),r(d+1:2*d),[0;0;1])];
    end

end
