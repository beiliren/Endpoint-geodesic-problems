% LSO3_demo2.m
clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The special orthogonal group SO(3): X^T*X = I_3, det(X) = 1
% Endow with left-invariant metric: <u,v> = tr(u^T*A*v), u,v\in the Lie
% algebra so(3), A is a fixed positive definite matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Auxiliary functions
options.A = diag([1,2,4]);
options.symm = @(M) (M+M')/2;
options.skew = @(M) (M-M')/2;
options.M_V = @(om) [0,-om(3),om(2);om(3),0,-om(1);-om(2),om(1),0];
options.iM_V = @(OM) [OM(3,2);OM(1,3);OM(2,1)];
% Bilinear form B_F(X1,X2) = A^-1*(A*X1\times X2)
options.B_F = @(X1,X2) options.A\cross(options.A*X1,X2);
options.TPj = @(X,Y) logm(X'*Y);

% initial conditions
X0 = eye(3);
angle_x = 80;
angle_y = 45;
XT = roty(angle_y)*rotx(angle_x);

T = 1;

options.eps0 = 10^-4;  % if two points are too close with respect to eps0
options.eps1 = 10^-4;  % iteration stop criteria
options.N = 5000;      % maximal iteration
step_s = [0.01,0.006,0.002,0.002,0.002,0.001];  % step size
n = [4,5,6,7,8,9];

% compare different algorithms
results_time = zeros(length(n),3);
for j = 1:length(n)
    ti = 0:1/n(j):1;
    X = [];
    for i = 1:n(j)+1
        X = [X,roty(ti(i)*angle_y)*rotx(ti(i)*angle_x)];
    end

    %% Riemannian gradient descent
    [~,tim_g,~] = Geodesic_LSO3_gradient(X,T,step_s(j),options);

    %% Leap frog
    [~,tim_lf,~] = Geodesic_LSO3_leapfrog(X,T,options);

    %% Newton's method
    [~,tim_n,~] = Geodesic_LSO3_newton(X,T,options);
    
    results_time(j,:) = [tim_g(end),tim_lf(end),tim_n(end)];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print results
fprintf('No.,\t RGD,\t Leapfrog,\t Newton\n');
for i = 1:length(n)
    fprintf('%d, %.4f, %.4f, %.4f\n', n(i), results_time(i,1),results_time(i,2),results_time(i,3));
end
