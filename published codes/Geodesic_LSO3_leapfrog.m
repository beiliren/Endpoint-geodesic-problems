% Geodesic_LSO3_leapfrog.m
function [Xnew,time_f,cost_f] = Geodesic_LSO3_leapfrog(X,T,options)


t_start = tic;
Xnew = X;
time_f = [];
cost_f = [];
[d,nd] = size(X);
n2 = nd/d;
n = n2-1; % number of segments
for k = 1:options.N
    V = zeros(2*d,n);
    for i = 1:n
        p = Xnew(:,d*(i-1)+1:d*i);
        q = Xnew(:,d*i+1:d*(i+1));
        sol1i = log_xy_LSO3(p,q,T/n,options);
        V(1:d,i) = sol1i.y(d*d+1:d*d+d,1); % v_i^+
        V(d+1:2*d,i) = sol1i.y(d*d+1:d*d+d,end); % v_{i+1}^-
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
    
    for i = 2:n
        yim1 = Xnew(:,d*(i-2)+1:d*(i-1));
        yip1 = Xnew(:,d*i+1:d*(i+1));
        sol1i = log_xy_LSO3(yim1,yip1,2*T/n,options);
        Xnew(:,d*(i-1)+1:d*i) = reshape(deval(sol1i,T/n,1:d*d),d,d);
    end

end
end
