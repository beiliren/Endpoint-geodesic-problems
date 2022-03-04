% Geodesic_LSO3_newton.m
function [Xnew, time_f, cost_f] = Geodesic_LSO3_newton(X,T,options)
% Input: 
% Output: 

t_start = tic;
Xnew = X;
time_f = [];
cost_f = [];
[d,nd] = size(X);
n1 = nd/d;
n = n1-1; % number of segments

for k = 1:options.N
    V = zeros(2*d,n);
    grad_V = zeros(4*d,d*n);
    for i = 1:n
        p = Xnew(:,d*(i-1)+1:d*i);
        q = Xnew(:,d*i+1:d*(i+1));
        sol1i = log_xy_LSO3(p,q,T/n,options);
        [dxlog1,dylog1,dxnlog1,dynlog1] = dlog_LSO3(p,q,sol1i,T/n,options);
        V(1:d,i) = sol1i.y(d*d+1:d*d+d,1);
        V(d+1:2*d,i) = sol1i.y(d*d+1:d*d+d,end);
        grad_V(1:d,d*(i-1)+1:d*i) = dxlog1; % d_{y_i} v_i^+
        grad_V(d+1:2*d,d*(i-1)+1:d*i) = dylog1; % d_{y_{i+1}} v_i^+
        grad_V(2*d+1:3*d,d*(i-1)+1:d*i) = -dxnlog1; % d_{y_{i+1}} -v_{i+1}^-
        grad_V(3*d+1:4*d,d*(i-1)+1:d*i) = -dynlog1; % d_{y_i} -v_{i+1}^-
    end
    
    F = V(1:d,2:end)-V(d+1:2*d,1:end-1);
    cur_er = sum(diag(F'*options.A*F))/2;
    time_f = [time_f,toc(t_start)];
    cost_f = [cost_f,cur_er];
    fprintf('%d, %.4f, %.4f\n',k, time_f(end), cost_f(end));
    if cur_er < options.eps1
        break;
    end
    
    vF = F(:);
    if norm(vF) > options.eps1
        dF = zeros(d*(n-1),d*(n-1));
        for j = 1:n-1
            if j == 1
                dF(d*(j-1)+1:d*j,d*(j-1)+1:d*j) = grad_V(1:d,d*j+1:d*(j+1))+grad_V(2*d+1:3*d,d*(j-1)+1:d*j);
                dF(d*(j-1)+1:d*j,d*j+1:d*(j+1)) = grad_V(d+1:2*d,d*j+1:d*(j+1));
            else
                if j < n-1
                    dF(d*(j-1)+1:d*j,d*(j-2)+1:d*(j-1)) = grad_V(3*d+1:4*d,d*(j-1)+1:d*j);
                    dF(d*(j-1)+1:d*j,d*(j-1)+1:d*j) = grad_V(1:d,d*j+1:d*(j+1))+grad_V(2*d+1:3*d,d*(j-1)+1:d*j);
                    dF(d*(j-1)+1:d*j,d*j+1:d*(j+1)) = grad_V(d+1:2*d,d*j+1:d*(j+1));
                else
                    dF(d*(j-1)+1:d*j,d*(j-2)+1:d*(j-1)) = grad_V(3*d+1:4*d,d*(j-1)+1:d*j);
                    dF(d*(j-1)+1:d*j,d*(j-1)+1:d*j) = grad_V(1:d,d*j+1:d*(j+1))+grad_V(2*d+1:3*d,d*(j-1)+1:d*j);
                end
            end
        end
        if rank(dF) == d*(n-1)
            v = linsolve(dF,-vF);
            for i1 = 1:n-1
                yi = Xnew(:,d*i1+1:d*(i1+1));
                Xnew(:,d*i1+1:d*(i1+1)) = exp_pv_LSO3(yi,v(d*(i1-1)+1:d*i1),options);
%                 yi = yi+(eye(d)-yi*yi'*options.A)*v(d*(i1-1)+1:d*i1);
%                 Xnew(:,i1+1) = yi/sqrt(yi'*options.A*yi);
            end
        else
            disp('Singular!');
        end
    else
        break;
    end

end
end