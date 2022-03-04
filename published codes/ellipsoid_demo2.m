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
n = [4,5,6,7,8,9];              % number of junctions: n-1
options.eps0 = 10^-4;  % if two points are too close with respect to eps0
options.eps1 = 10^-4;  % iteration stop criteria
options.N = 5000;      % maximal iteration
step_s = [0.01,0.001,0.001,0.001,0.001,0.001];     % step size


% compare different algorithms
results_time = zeros(length(n),3);
for j = 1:length(n)
    thetai = theta0:(theta1-theta0)/n(j):theta1;
    phii = phi0:(phi1-phi0)/n(j):phi1;
    X = [];
    for i = 1:n(j)+1
        X = [X,ellipsoid_s(thetai(i),phii(i))];
    end

    %% Riemannian gradient descent
    [~,tim_g,~,cost_f_g] = Geodesic_ellipsoid_gradient(X,T,step_s(j),options);

    %% Leap frog
    [~,tim_lf,cost_f_lf] = Geodesic_ellipsoid_leapfrog(X,T,options);

    %% Newton's method
    [~,tim_n,cost_f_n] = Geodesic_ellipsoid_newton(X,T,options);
    
    results_time(j,:) = [tim_g(end),tim_lf(end),tim_n(end)];

    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print results
fprintf('RGD,\t Leapfrog,\t Newton\n');
for i = 1:length(n)
    fprintf('%d, %.4f, %.4f, %.4f\n', n(i), results_time(i,1),results_time(i,2),results_time(i,3));
end


