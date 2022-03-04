% log_xy_LSO3.m
function sol = log_xy_LSO3(x,y,T,options)
% For simplicity, we denote rotation matrix by its Euler angles.
% Assume x = Rz(t3)*Ry(t2)*Rx(t1)
%          =[cos(t2)*cos(t3) -cos(t1)*sin(t3)+sin(t1)*sin(t2)*cos(t3) sin(t1)*sin(t3)+cos(t1)*sin(t2)*cos(t3)
%            cos(t2)*sin(t3) cos(t1)*cos(t3)+sin(t1)*sin(t2)*sin(t3) -sin(t1)*cos(t3)+cos(t1)*sin(t2)*sin(t3)
%            -sin(t2) sin(t1)*cos(t2) cos(t1)*cos(t2)]
%    dx = x*hv -->  
%    d theta = [1 sin(theta1)*tan(theta2) cos(theta1)*tan(theat2)
%               0 cos(theta1) -sin(theta1)
%               0 sin(theta1)/cos(theta2) cos(theta1)/cos(theta2)]*v

eps0 = options.eps0;
A = options.A;
iM_V = options.iM_V;
M_V = options.M_V;
B_F = options.B_F;
TPj = options.TPj;

d = size(x,1);
if norm(x-y) < eps0
%     sol.y = [x(:);iM_V(TPj(x,y)./T)]; 
%     sol.y = [x(1:d*(d-1));iM_V(TPj(x,y)./T)]; 
    sol.y = [x(:);reshape(TPj(x,y)./T,[],1)]; 
else
    tmesh = linspace(0,T,21);
    solinit = bvpinit(tmesh,@ini_guess);
    bvpfcn = @(t,y) geodesiceq_LSO3(t,y);
    opts = bvpset('FJacobian',@jac);
    sol = bvp4c(bvpfcn,@bcfcn,solinit,opts); % bvp5c
    for i = 1:length(sol.x)
        skmv = reshape(sol.y(d*d+1:2*d*d,i),d,d);
        skmv = (skmv-skmv')./2;
        sol.y(d*d+1:d*d+d,i) = iM_V(skmv);
    end

end

    % initial guess
    function g = ini_guess(tt)
%         g = [x(:);tt*iM_V(TPj(x,y)./T)];
%         g = [x(1:d*(d-1));tt*iM_V(TPj(x,y)./T)];
        g = [x(:);reshape(tt.*TPj(x,y)./T,d*d,1)];
    end

    % geodesic equation
    function drdt = geodesiceq_LSO3(t,r)
%         Rv = reshape(r,d,[]);
%         R = Rv(:,1:d);
%         v = Rv(:,d+1);
%         dRdt = R*M_V(v);
%         drdt = [dRdt(:);B_F(v,v)];
        
%         R12v = reshape(r,d,d);
%         R1 = R12v(:,1);
%         R2 = R12v(:,2);
%         R = [R1,R2,cross(R1,R2)];
%         v = R12v(:,d);
%         dR12dt = R*M_V(v);
%         drdt = [reshape(dR12dt(:,1:2),d*(d-1),1);B_F(v,v)];
        
        RV = reshape(r,d,2*d);
        R = RV(:,1:d);
        V = RV(:,d+1:2*d);
        V = (V-V')./2;
        v = iM_V(V);
        dRdt = R*M_V(v);
        drdt = [dRdt(:);reshape(M_V(B_F(v,v)),d*d,1)];
    end

    % boundary conditions
    function res = bcfcn(r0,r1)
%         res = [r0(1:d*(d-1))-reshape(x(:,1:d-1),d*(d-1),1)
%             r1(1:d*(d-1))-reshape(y(:,1:d-1),d*(d-1),1)];
        
        res = [r0(1:d*d)-x(:)
            r1(1:d*d)-y(:)];
    end

   % Jacobian
    function dydt = jac(t,y)
        % y is 12*1 vector
%         yR = reshape(y(1:d*d),d,d);
%         yv = y(d*d+1:d*d+d);
%         dRdv = [-yR*M_V([1,0,0]);-yR*M_V([0,1,0]);-yR*M_V([0,0,1])];
%         dydt = [-kron(M_V(yv),eye(d)),dRdv;zeros(d,d*d),M_V(A*yv)-M_V(yv)*A];
        % y is 18*1 vector
%         yR = reshape(y(1:d*d),d,d);
%         yV = reshape(y(d*d+1:2*d*d),d,d);
%         yV = (yV-yV')./2;
%         yv = iM_V(yV);
%         dydt = [-kron(M_V(yv),eye(d)),kron(eye(d),yR);zeros(d*d,d*d),kron(eye(d),M_V(yv)-A\M_V(yv)*A)-kron(M_V(yv),eye(d))+kron(M_V(yv)*A,inv(A))];

        yR = reshape(y(1:d*d),d,d);
        yV = reshape(y(d*d+1:2*d*d),d,d);
        yV = (yV-yV')./2;
        dRVdR = kron(yV',eye(d));
        dRVdV = kron(eye(d),yR);
        dfVdV = dRVdR+kron(eye(d),yV-A\yV*A)-kron(yV'*A,inv(A));
        dydt = [dRVdR,dRVdV;zeros(d*d,d*d),dfVdV];
    end
end


function theta = x_theta(x)
% map x as theta
theta2 = asin(-x(3,1));
if x(1,1) ~= 0
    theta3 = atan(x(2,1)/x(1,1));
    theta1 = asin(x(3,2)/cos(theta2));
else
    if cos(theta2) ~= 0
        theta3 = pi/2;
        theta1 = asin(x(3,2)/cos(theta2));
    else
        fprintf('Have trouble here! theta1-theta3=%.4f\n',asin(x(1,2)));
    end
end
theta = [theta1;theta2;theta3];
end

function x = theta_x(theta)
% map theta as x
x = [cos(theta(2))*cos(theta(3)) -cos(theta(1))*sin(theta(3))+sin(theta(1))*sin(theta(2))*cos(theta(3)) sin(theta(1))*sin(theta(3))+cos(theta(1))*sin(theta(2))*cos(theta(3))
     cos(theta(2))*sin(theta(3)) cos(theta(1))*cos(theta(3))+sin(theta(1))*sin(theta(2))*sin(theta(3)) -sin(theta(1))*cos(theta(3))+cos(theta(1))*sin(theta(2))*sin(theta(3))
     -sin(theta(2)) sin(theta(1))*cos(theta(2)) cos(theta(1))*cos(theta(2))];
end