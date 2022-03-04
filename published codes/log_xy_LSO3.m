% log_xy_LSO3.m
function sol = log_xy_LSO3(x,y,T,options)

eps0 = options.eps0;
A = options.A;
iM_V = options.iM_V;
M_V = options.M_V;
B_F = options.B_F;
TPj = options.TPj;

d = size(x,1);
if norm(x-y) < eps0
    sol.y = [x(:);reshape(TPj(x,y)./T,[],1)]; 
else
    tmesh = linspace(0,T,21);
    solinit = bvpinit(tmesh,@ini_guess);
    bvpfcn = @(t,y) geodesiceq_LSO3(t,y);
    opts = bvpset('FJacobian',@jac);
    sol = bvp4c(bvpfcn,@bcfcn,solinit,opts); % bvp5c
    for i = 1:length(sol.x)
        skmv = reshape(sol.y(d*d+1:2*d*d,i),d,d);
        skmv = options.skew(skmv);
        sol.y(d*d+1:d*d+d,i) = iM_V(skmv);
    end

end

    % initial guess
    function g = ini_guess(tt)
        g = [x(:);reshape(tt.*TPj(x,y)./T,d*d,1)];
    end

    % geodesic equation
    function drdt = geodesiceq_LSO3(t,r)
        RV = reshape(r,d,2*d);
        R = RV(:,1:d);
        V = RV(:,d+1:2*d);
        V = options.skew(V);
        v = iM_V(V);
        dRdt = R*M_V(v);
        drdt = [dRdt(:);reshape(M_V(B_F(v,v)),d*d,1)];
    end

    % boundary conditions
    function res = bcfcn(r0,r1)
        res = [r0(1:d*d)-x(:)
               r1(1:d*d)-y(:)];
    end

   % Jacobian
    function dydt = jac(t,y) 
        yR = reshape(y(1:d*d),d,d);
        yV = reshape(y(d*d+1:2*d*d),d,d);
        yV = (yV-yV')./2;
        dRVdR = kron(yV',eye(d));
        dRVdV = kron(eye(d),yR);
        dfVdV = dRVdR+kron(eye(d),yV-A\yV*A)-kron(yV'*A,inv(A));
        dydt = [dRVdR,dRVdV;zeros(d*d,d*d),dfVdV];
    end
end
