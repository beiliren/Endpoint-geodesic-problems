%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Methods for solving endpoint geodesic problems proposed in
% 
% Lyle Noakes and Erchuan Zhang, 
%        "Finding geodesics joining given points"
% Advances in Computational Mathematics, 2022.
% 
% Copyright (c) 2022 Erchuan Zhang (erchuan.zhang@ecu.edu.au)
% School of Science, Edith Cowan University, Australia
% 
% Please acknowledge the authors by citing the above paper in any academic 
% publications that have made use of this package or part of it.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ellipsoid: x^2/a^2+y^2/b^2+z^2/c^2=1 
% parametrization: x = a*cos(phi)*sin(theta)
%                  y = b*sin(phi)*sin(theta)
%                  z = c*cos(theta)  phi\in [0,2*pi), theta\in [0,pi]
% Riemannian metric on ellipsoid is induced from Euclidean metric
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameterization of the 2-ellipsoid in R^3
a = 6; b = 2; c = 1;
A = diag([1/a^2,1/b^2,1/c^2]);
options.A = A;
ellipsoid_s = @(theta,phi) [a*cos(phi)*sin(theta);b*sin(phi)*sin(theta);c*cos(theta)];

% Christoffel symbol C_S(x,u,v)=\Gamma_x(u,v)
options.C_S = @(x,u,v) ((u'*A*v)/(x'*A^2*x))*A*x;

% derivative of Christoffel symbol DC_S(x,w,u)=d\Gamma_x(w)(u,u)
options.DC_S = @(x,w,u) ((u'*A*u)/(x'*A^2*x)).*A*w-2*((u'*A*u*x'*A^2*w)/(x'*A^2*x)^2).*A*x;

% projection vector y-x as a tangent vector at x
options.TPj = @(x,y,T) (y-x-((y-x)'*A*x)/(x'*A^2*x)*A*x)/T;

% exponential map
options.retr = @(x,v,ops) exp_pv_ellipsoid(x,v,ops);

% initial conditions
theta0 = pi/2;
phi0 = 0;
theta1 = pi/6;
phi1 = 2*pi/3;
x0 = ellipsoid_s(theta0,phi0);  % initial position
x1 = ellipsoid_s(theta1,phi1);  % final position
T = 1;                          % time interval
n = 5;                          % number of midpoints: n-1
options.eps0 = 10^-4;  % if two points are too close with respect to eps0
options.eps1 = 10^-4;  % iteration stop criteria
options.N = 5000;      % maximal iteration


% Solve the boundary value problem directly
% sol10 = log_xy_ellipsoid(x0,x1,T,options);
% fprintf('The initial point is: [%.4f, %.4f, %.4f]\n', x0');
% fprintf('The final point is: [%.4f, %.4f, %.4f]\n', x1');
% fprintf('The final point calculated by shooting is: [%.4f, %.4f, %.4f]\n', sol10.y(1:3,end)');

% initial midpoints
thetai = theta0:(theta1-theta0)/n:theta1;
phii = phi0:(phi1-phi0)/n:phi1;
X = [];
for i = 1:n+1
    X = [X,ellipsoid_s(thetai(i),phii(i))];
end

% update midpoints
%% Newton's method
[Xnew,tim,cost_f] = Geodesic_ellipsoid_newton(X,T,options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot results
figure;
% plot initial midpoints
plot3(X(1,2:n),X(2,2:n),X(3,2:n),'mo','MarkerFaceColor','m');
hold on;

% plot initial piecewise geodesic
geo_initial = [];
for i = 1:n
    solge = log_xy_ellipsoid(X(:,i),X(:,i+1),T/n,options);
    geo_initial = [geo_initial,solge.y(1:3,2:end-1)];
end
plot3(geo_initial(1,:),geo_initial(2,:),geo_initial(3,:),'g-','LineWidth',2);

% plot final midpoints
plot3(Xnew(1,2:n),Xnew(2,2:n),Xnew(3,2:n),'bo','MarkerFaceColor','b');

% plot final generated geodesic
geo_final = [];
for i = 1:n
    solge = log_xy_ellipsoid(Xnew(:,i),Xnew(:,i+1),T/n,options);
    geo_final = [geo_final,solge.y(1:3,2:end-1)];
end
plot3(geo_final(1,:),geo_final(2,:),geo_final(3,:),'r-','LineWidth',2);

% plot two given points
plot3(X(1,1),X(2,1),X(3,1),'ko',X(1,n+1),X(2,n+1),X(3,n+1),'ko','MarkerFaceColor','k');
text(X(1,1),X(2,1),X(3,1)-0.1,'$x_0$','Interpreter','latex');
text(X(1,n+1),X(2,n+1),X(3,n+1)+0.1,'$x_1$','Interpreter','latex');

% plot ellipsoid
[xa,xb,xc] = ellipsoid(0,0,0,a,b,c);
lightGrey = 0.8*[1,1,1];
surface(xa,xb,xc,'FaceColor','none','EdgeColor',lightGrey);
set(gca,'visible','off')
hold off;
