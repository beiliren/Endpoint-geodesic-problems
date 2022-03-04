% dlog_ellipsoid.m
function [dxlog,dylog,dxnlog,dynlog] = dlog_ellipsoid(x,y,solin,T,options)
% Input: x,y: two points on M (vectors of d*1)
%        C_S: Christoffel symbol on M, C_S(x,v,vv)=\Gamma_x(v,vv)
%        DC_S: derivative of Christoffel symbol, DC_S(x,J,v)=d\Gamma_x(J)(v,v)
%        eps0: small positive value
% Output: dxlog: d_xlog(x,y) along the direction e_i
%         dylog: d_ylog(x,y) along the direction e_i

A = options.A;
eps0 = options.eps0;

d = size(x,1); % dimension of M
if sqrt((x-y)'*A*(x-y)) < eps0 %norm(x-y) < eps0
    dxlog = -eye(d);
    dylog = eye(d);
    dxnlog = -eye(d);
    dynlog = eye(d);
else
    e0 = zeros(2*d,1);
    dxexp = zeros(d,d);
    dJ1 = zeros(d,d);
    for i = 1:d
        e0(i,1) = 1;
        [~,sol1] = ode45(@jacobieq,[0,T],e0);
        dxexp(:,i) = sol1(end,1:d)'; % J_1(1) with J_1(0)=e,dJ_1(0)=0
        dJ1(:,i) = sol1(end,d+1:2*d)';
        e0(i,1) = 0;
    end
    dvexp = zeros(d,d);
    dJ2 = zeros(d,d);
    for i = d+1:2*d
        e0(i,1) = 1;
        [~,sol2] = ode45(@jacobieq,[0,T],e0);
        dvexp(:,i-d) = sol2(end,1:d)'; % J_2(1) with J_2(0)=0,dJ_2(0)=e
        dJ2(:,i-d) = sol2(end,d+1:2*d)';
        e0(i,1) = 0;
    end
    dylog = inv(dvexp);
    dxlog = -dylog*dxexp;
    M = [dxexp,dvexp;dJ1,dJ2];
    invM = inv(M);
    dynlog = inv(invM(1:d,d+1:2*d));
    dxnlog = -dynlog*invM(1:d,1:d);
end

    function dJdt = jacobieq(t,J)
        % create the Jacobi equation
        rr = deval(solin,t,1:d);
        drr = deval(solin,t,d+1:2*d);
        dJdt =[J(d+1:2*d);
               -2.*options.C_S(rr,drr,J(d+1:2*d))-options.DC_S(rr,J(1:d),drr)];
    end

end