% exp_pv_LSO3.m
function y = exp_pv_LSO3(x,v,options)

d = size(x,1);
if norm(v) <= options.eps0
    y = x;
else
    odefun = @(t,r) geodesiceq_G(t,r);
    v0 = [x,v];
    [~,sol] = ode45(odefun,[0,1],v0(:));
    rr1 = sol(end,1:d*d);
    y = reshape(rr1,[d,d]);
end

% geodesiceq_G.m
   function drdt = geodesiceq_G(t,r)
       RV = reshape(r,d,d+1);
       w = RV(:,d+1);
       RRV = [RV(:,1:d)*options.M_V(w),options.B_F(w,w)];
       drdt = RRV(:);
   end

end