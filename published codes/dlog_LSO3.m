% dlog_LSO3.m
function [dxlog,dylog,dxnlog,dynlog] = dlog_LSO3(x,y,solin,T,options)

eps0 = options.eps0;
B_F = options.B_F;

d = size(x,1); % dimension of M
if norm(x-y) <= eps0
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
        dxexp(:,i) = sol1(end,1:d)';
        dJ1(:,i) = sol1(end,d+1:2*d)';
        e0(i,1) = 0;
    end
    dvexp = zeros(d,d);
    dJ2 = zeros(d,d);
    for i = d+1:2*d
        e0(i,1) = 1;
        [~,sol2] = ode45(@jacobieq,[0,T],e0);
        dvexp(:,i-d) = sol2(end,1:d)';
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

    function drJdt = jacobieq(t,rJ)
        v = deval(solin,t,d*d+1:d*d+d);
        drJdt = [rJ(d+1:2*d);
            cross(rJ(d+1:2*d),v)+...
            B_F(rJ(d+1:2*d),v)+...
            B_F(v,rJ(d+1:2*d))+...
            cross(rJ(1:d),B_F(v,v))- ...
            B_F(cross(rJ(1:d),v),v)- ...
            B_F(v,cross(rJ(1:d),v))];
    end

end