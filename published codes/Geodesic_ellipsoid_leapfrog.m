% Geodesic_ellipsoid_leapfrog.m
function [Xnew,tim,cost_f] = Geodesic_ellipsoid_leapfrog(X,T,options)
% Input: X: initial midpoints on M (d*(n+1), first column is x0, last column is x1)
%        N: maximal number of iterations
%        C_S: Christoffel symbol
%        DC_S: derivative of Christoffel symbol, DC_S(x,J,v)=d\Gamma_x(J)(v,v)
%        TPj: projection from the embedded Euclidean space to TM
%        step_s: step size for iteration
%        eps0: small positive value
%        eps1: stop criteria for iteration

t_start = tic;
Xnew = X(:,2:end-1);
tim = [];
cost_f = [];
[d,n_1] = size(Xnew); % n_1: number of junctions
n = n_1+1; % n: number of segments
for k = 1:options.N
    cur_er = 0;
    for i = 1:n_1
        if i == 1
            yim1 = X(:,1);
        else
            yim1 = Xnew(:,i-1);
        end
        if i < n_1
            yip1 = Xnew(:,i+1);
        else
            yip1 = X(:,end);
        end
        yi = Xnew(:,i);
        sol1i = log_xy_ellipsoid(yim1,yi,T/n,options);
        sol2i = log_xy_ellipsoid(yi,yip1,T/n,options);
        vip = sol2i.y(d+1:2*d,1);
        vim = sol1i.y(d+1:2*d,end);
        cur_er = cur_er+norm(vip-vim)^2;
    end
    tim = [tim,toc(t_start)];
    cost_f = [cost_f,cur_er/2];
    fprintf('%d, %.4f, %.4f\n',k, tim(end), cost_f(end));
    if cur_er/2 < options.eps1 
        Xnew = [X(:,1),Xnew,X(:,end)];
        break;
    end
    
    
    for i = 1:n_1
        if i == 1
            yim1 = X(:,1);
        else
            yim1 = Xnew(:,i-1);
        end
        if i < n_1
            yip1 = Xnew(:,i+1);
        else
            yip1 = X(:,end);
        end
        sol1i = log_xy_ellipsoid(yim1,yip1,2*T/n,options);
        Xnew(:,i) = deval(sol1i,T/n,1:d);
        Xnew(:,i) = Xnew(:,i)/sqrt(Xnew(:,i)'*options.A*Xnew(:,i));
    end

end
end
