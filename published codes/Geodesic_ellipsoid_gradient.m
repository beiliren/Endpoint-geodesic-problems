% Geodesic_ellipsoid_gradient.m
function [Xnew,tim,norm_df,cost_f] = Geodesic_ellipsoid_gradient(X,T,step_s,options)
% Input: X: initial midpoints on M (d*(n+1), first column is x0, last column is x1)
%        N: maximal number of iterations
%        C_S: Christoffel symbol
%        DC_S: derivative of Christoffel symbol, DC_S(x,J,v)=d\Gamma_x(J)(v,v)
%        TPj: projection from the embedded Euclidean space to TM
%        step_s: step size for iteration
%        eps0: small positive value
%        eps1: stop criteria for iteration

t_start = tic;
Xnew = X;
tim = [];
norm_df = [];
cost_f = [];
[d,nd] = size(X);
n = nd-1; % number of segments
for k = 1:options.N
    V = zeros(2*d,n);
    grad_V = zeros(4*d,d*n);
    for i = 1:n
        p = Xnew(:,i);
        q = Xnew(:,i+1);
        sol1i = log_xy_ellipsoid(p,q,T/n,options);
        [dxlog1,dylog1,dxnlog1,dynlog1] = dlog_ellipsoid(p,q,sol1i,T/n,options);
        V(1:d,i) = sol1i.y(d+1:2*d,1); % v_i^+
        V(d+1:2*d,i) = sol1i.y(d+1:2*d,end); % v_{i+1}^-
        grad_V(1:d,d*(i-1)+1:d*i) = dxlog1; % d_{y_i} v_i^+
        grad_V(d+1:2*d,d*(i-1)+1:d*i) = dylog1; % d_{y_{i+1}} v_i^+
        grad_V(2*d+1:3*d,d*(i-1)+1:d*i) = -dxnlog1; % d_{y_{i+1}} -v_{i+1}^-
        grad_V(3*d+1:4*d,d*(i-1)+1:d*i) = -dynlog1; % d_{y_i} -v_{i+1}^-
    end
    c_f = 1/2*sum(vecnorm(V(1:d,2:end)-V(d+1:2*d,1:end-1)).^2);
    tim = [tim,toc(t_start)];
    cost_f = [cost_f,c_f];
    fprintf('%d, %.4f, %.4f\n',k, tim(end), cost_f(end));
    if c_f < options.eps1
        break;
    end
    
    % add time limitation
    if tim(end) > 100
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
        % Riemannian gradient
        p = Xnew(:,j+1);
        df(:,j) = df(:,j)-p*p'*options.A*df(:,j);
%         Xnew(:,j+1) = exp_pv_ellipsoid(p,-step_s*df(:,j),options);
        Xnew(:,j+1) = (p-step_s*df(:,j))./sqrt((p-step_s*df(:,j))'*options.A*(p-step_s*df(:,j)));
    end
    n_df = sum(vecnorm(df))/(n-1);
    norm_df = [norm_df,n_df];

end
end
