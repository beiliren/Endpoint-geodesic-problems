% LSO3_demo1.m
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
% guess for initial velocity
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

% initial midpoints
n = 5;        % number of midpoints: n-1
ti = 0:1/n:1;
Xini = [];
for i = 1:n+1
    Xini = [Xini,roty(ti(i)*angle_y)*rotx(ti(i)*angle_x)];
end


% update midpoints
%% Newton's method
[Xnew,time_f,cost_f] = Geodesic_LSO3_newton(Xini,T,options); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot sphere
figure;
[a,b,c] = sphere(20);
lightGrey = 0.8*[1,1,1];
surface(a,b,c,'FaceColor','none','EdgeColor',lightGrey);
set(gca,'visible','off')
hold on;

% plot initial rotation matrix
X00 = [X0(:,1),zeros(3,1),X0(:,2),zeros(3,1),X0(:,3)];
plot3(X00(1,:),X00(2,:),X00(3,:),'r-','LineWidth',2); %hold on;
% plot final rotation matrix
X10 = [XT(:,1),zeros(3,1),XT(:,2),zeros(3,1),XT(:,3)];
plot3(X10(1,:),X10(2,:),X10(3,:),'k-','LineWidth',2); %hold off;

% plot trajectory of initial midpoints
plot3(Xini(1,:),Xini(2,:),Xini(3,:),'mo','MarkerFaceColor','m');

% plot final midpoints
plot3(Xnew(1,:),Xnew(2,:),Xnew(3,:),'bo','MarkerFaceColor','b');

% plot final generated geodesic
geo_final = [];
for i = 1:n
    ge = log_xy_LSO3(Xnew(:,3*(i-1)+1:3*i),Xnew(:,3*i+1:3*(i+1)),T/n,options);
    geo_final = [geo_final,reshape(ge.y(1:9,2:end-1),3,[])];
end
plot3(geo_final(1,:),geo_final(2,:),geo_final(3,:),'g.','MarkerSize',10);
hold off;
