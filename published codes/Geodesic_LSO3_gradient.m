% Geodesic_LSO3_gradient.m
function [Xnew,time_f,cost_f] = Geodesic_LSO3_gradient(X,T,step_s,options)

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
    
    dv = V(1:d,2:end)-V(d+1:2*d,1:end-1);
    cur_er = sum(diag(dv'*options.A*dv))/2;
    time_f = [time_f,toc(t_start)];
    cost_f = [cost_f,cur_er];
    fprintf('%d, %.4f, %.4f\n',k, time_f(end), cost_f(end));
    if cur_er < options.eps1
        break;
    end
    
    % add time limitation
    if time_f(end) > 100
        break;
    end
    
    df = zeros(d,n-1); % gradient of cost function
    for j = 1:n-1
        df(:,j) = (grad_V(1:d,d*j+1:d*(j+1))+grad_V(2*d+1:3*d,d*(j-1)+1:d*j))'*(V(1:d,j+1)-V(d+1:2*d,j));
        % % A*
        if n > 2
            if j == 1
                df(:,j) = df(:,j)+grad_V(3*d+1:4*d,d*j+1:d*(j+1))'*(V(1:d,j+2)-V(d+1:2*d,j+1));
            else
                if j < n-1
                    df(:,j) = df(:,j)+grad_V(3*d+1:4*d,d*j+1:d*(j+1))'*(V(1:d,j+2)-V(d+1:2*d,j+1));
                    df(:,j) = df(:,j)+grad_V(d+1:2*d,d*(j-1)+1:d*j)'*(V(1:d,j)-V(d+1:2*d,j-1));
                else
                    df(:,j) = df(:,j)+grad_V(d+1:2*d,d*(j-1)+1:d*j)'*(V(1:d,j)-V(d+1:2*d,j-1));
                end
            end
        end
        p = Xnew(:,d*j+1:d*(j+1));
        Xnew(:,d*j+1:d*(j+1)) = exp_pv_LSO3(p,-step_s*df(:,j),options);
    end

end
end
