close all
tic
%================= Choice of nodes ====================
%%1. Halton nodes
% load('100HaltonNodesb2b3.mat'); 
% nodes=[nodes(:,1), nodes(:,2), nodes(:,3)];
% n=length(nodes);

%%2. Spiral nodes
load('500spiralpoints.mat');
nodes=[points(:,1), points(:,2), points(:,3)];
n=length(nodes);

%%Delaunay triangulation of the sphere
[face_num, face] = sphere_delaunay ( n, nodes');
triangles = face';

%===============Choice of evaluation points================
%%1. Spiral points
% load('600spiralpoints.mat');
% x = [points(:,1), points(:,2), points(:,3)];
% points=x;
% nr_points = length(points);

%%2. Random points
load('1000points-random.mat')
random_points=points;
TH_points = 2*pi*random_points;
PH_points = asin(-1+2*random_points);
[xP,yP,zP] = sph2cart(TH_points,PH_points,1);
points = [xP' yP' zP'];
nr_points = length(points);
%==============================================

miu=2;
results_for_rmse_1= zeros(1, nr_points);
results_for_mae_1 = zeros(1, nr_points);
results_for_rmse_2= zeros(1, nr_points);
results_for_mae_2 = zeros(1, nr_points);
true_vals = zeros(1, nr_points);

for j=1:nr_points
    a1=0; b=0; a2=0;
    [a1,a2]=shepard_bern(points(j,:), n, nodes, triangles, miu);
    b=fc_3d(points(j,:));
    results_for_mae_1(j) = abs(a1-b);
    results_for_mae_2(j) = abs(a2-b);
    results_for_rmse_1(j)=a1;
    results_for_rmse_2(j)=a2;
    true_vals(j)=b;
end
%====================================================

max_abs_err_order_1 = max(results_for_mae_1)
mean_abs_err_order_1 = mean(results_for_mae_1)
rmse_order_1=(sum((results_for_rmse_1-true_vals).^2)/nr_points)^1/2

max_abs_err_order_2 = max(results_for_mae_2)
mean_abs_err_order_2 = mean(results_for_mae_2)
rmse_order_2=(sum((results_for_rmse_2-true_vals).^2)/nr_points)^1/2

toc
%====================================================

function [s1,s2] = shepard_bern(x, n, nodes, triangles, miu)
m = length(triangles);
s1=0; s2=0;
 v = sum_all_triangles(x, nodes, triangles, miu);
    for i=1:m
        phi_i = phi(x, i ,nodes, triangles, miu);
    %=======================================================
        node1 = nodes(triangles(i,1),:); 
        node2 = nodes(triangles(i,2),:);
        node3 = nodes(triangles(i,3),:);
   %=======================================================
        spher_area = spherical_area(node1, node2,node3);
        lambda1 = spherical_area(x,node2,node3)/spher_area;
        lambda2 = spherical_area(node1,x,node3)/spher_area;
        lambda3 = spherical_area(node1,node2,x)/spher_area;
    %=======================================================
        xx = atan4(x(2),x(1));
        yy = arc_cosine(x(3)); 
        x1=atan4(node1(2),node1(1));
        y1=arc_cosine(node1(3));
        node1=[x1 y1];
        x2=atan4(node2(2),node2(1));
        y2=arc_cosine(node2(3));
        node2=[x2 y2];
        x3=atan4(node3(2),node3(1));
        y3=arc_cosine(node3(3));
        node3=[x3 y3];
        point = [xx yy];
        res1 = bern_ord_1(node1, node2, node3, lambda1, lambda2, lambda3);
        res2 = bern_ord_2(node1, node2, node3, lambda1, lambda2, lambda3);
    %=======================================================
        s1 = s1+phi_i*res1/v;
        s2 = s2+phi_i*res2/v;
    end
end

function tr_vertex_i = triangles_with_vertex_i(triangles,i)
    [row,~]=find(triangles==i);
    tr_vertex_i = triangles(row,:);
end

function f =  bern_ord_1(node1, node2, node3, lambda1, lambda2, lambda3)
    f = 0;
    f = lambda1*fc_2d(node1) +  lambda2*fc_2d(node2)+ lambda3*fc_2d(node3);
end

function f =  bern_ord_2(node1, node2, node3, lambda1, lambda2, lambda3)
    f = 0;
    x1=node1(1);
    x2=node2(1);
    x3=node3(1);
    y1=node1(2);
    y2=node2(2);
    y3=node3(2);
    f = fc_2d(node1)* lambda1 + fc_2d(node2)* lambda2 + fc_2d(node3)* lambda3...
    + 0.5*lambda1*lambda2 *((x1-x2)*(der_1x(node2)-der_1x(node1)) + (y1-y2)*(der_1y(node2)-der_1y(node1))) ...
    + 0.5*lambda1*lambda3 *((x3-x1)*(der_1x(node1)-der_1x(node3)) + (y3-y1)*(der_1y(node1)-der_1y(node3))) ...
    + 0.5*lambda2*lambda3 *((x2-x3)*(der_1x(node3)-der_1x(node2)) + (y2-y3)*(der_1y(node3)-der_1y(node2)));
end

function r = phi(x, i ,nodes, tr, miu)
    numer = 1/prod([geodesic_dist(x,nodes(tr(i,1),:)), geodesic_dist(x,nodes(tr(i,2),:)), geodesic_dist(x,nodes(tr(i,3),:))])^miu;
    r = numer;
end 

function s = sum_all_triangles(x, nodes, tr, miu)
    s = 0;
    m = length(tr);
    for i =1:m
        s = s+  1/prod([geodesic_dist(x,nodes(tr(i,1),:)), geodesic_dist(x,nodes(tr(i,2),:)), geodesic_dist(x,nodes(tr(i,3),:))])^miu;
    end
end

function g = geodesic_dist(x,y)
    g = arc_cosine(x*(y'));
end


function rez = spherical_area(a,b,c)
 s = a* cross(b, c)';
 s = s/( 1 + a*b' + b*c' + a*c');
 rez = 2*atan(s);
end

function f = fc_3d(u)          %%3D
x = u(1); y = u(2); z=u(3);

%===================== f1 =============================
f = 0.75*exp(-((9*x-2)^2)/4 - ((9*y-2)^2)/4 - ((9*z-2)^2)/4)  ...
+ 0.75*exp(-((9*x+1)^2)/49 - ((9*y+1)^2)/10 - ((9*z+1)^2)/10) ...
+ 0.5*exp(-((9*x-7)^2)/4 - ((9*y-3)^2)/4 - ((9*z-5)^2)/4) ...
- 0.2*exp (-((9*x-4))^2 - (9*y-7)^2 - (9*z-5)^2); %%f1 

%===================== f2 =============================
% f = (1.25+cos(5.4*y))*cos(6*z)/(6+6*(3*x-1)^2); %%f2

%===================== f3 =============================
% f = exp(-(81/16)*((x-0.5)^2+(y-0.5)^2+(z-0.5)^2))/3; %%f3

%===================== f4 =============================
% f =  exp(-(81/4)*((x-0.5)^2+(y-0.5)^2+(z-0.5)^2))/3; %%f4 

%===================== f5 =============================
% f = 0.1*(exp(x) + exp(y+z));                      %%f5

end

function f = fc_2d(u)            %%2D
phi = u(1); theta = u(2);
x = sin(theta)*cos(phi);
y = sin(theta)*sin(phi);
z = cos(theta);

%===================== f1 =============================
f = 0.75*exp(-((9*x-2)^2)/4 - ((9*y-2)^2)/4 - ((9*z-2)^2)/4)  ...
+ 0.75*exp(-((9*x+1)^2)/49 - ((9*y+1)^2)/10 - ((9*z+1)^2)/10) ...
+ 0.5*exp(-((9*x-7)^2)/4 - ((9*y-3)^2)/4 - ((9*z-5)^2)/4) ...
- 0.2*exp (-((9*x-4))^2 - (9*y-7)^2 - (9*z-5)^2); %%f1

%===================== f2 =============================
% f = (1.25+cos(5.4*y))*cos(6*z)/(6+6*(3*x-1)^2); %%f2 

% ===================== f3 =============================
% f = exp(-(81/16)*((x-0.5)^2+(y-0.5)^2+(z-0.5)^2))/3; %%f3 

%===================== f4 =============================
% f =  exp(-(81/4)*((x-0.5)^2+(y-0.5)^2+(z-0.5)^2))/3; %%f4 

%===================== f5 =============================
% f = 0.1*(exp(x) + exp(y+z));                      %%f5

end

function f = der_1x(u)
phi=u(1); theta = u(2);
x = sin(theta)*cos(phi);
y = sin(theta)*sin(phi);
z = cos(theta);

%===================== f1 =============================
A = 0.75*exp(-((9*x-2)^2)/4 - ((9*y-2)^2)/4 - ((9*z-2)^2)/4);
B =  0.75*exp(-((9*x+1)^2)/49 - ((9*y+1)^2)/10 - ((9*z+1)^2)/10);
C = 0.5*exp(-((9*x-7)^2)/4 - ((9*y-3)^2)/4 - ((9*z-5)^2)/4) ;
D = - 0.2*exp (-(9*x-4)^2 - (9*y-7)^2 - (9*z-5)^2);

fx = -0.5*A*(9*x-2) -2/49*B*(9*x+1) -0.5*C*(9*x-7) -2*D*(9*x-4);
fy = -0.5*A*(9*y-2) -0.2*B*(9*y+1) -0.5*C*(9*y-3) -2*D*(9*y-7);

f = -fx*sin(phi)*sin(theta) + fy*sin(theta)*cos(phi);  %%f1

% % %===================== f2 =============================
% A = (6+6*(3*x-1)^2)^(-1);
% B = (1.25 + cos(5.4*y));
% C = cos(6*z);
% 
% f=(6+6*(3*x-1)^2)^(-2)*(36*x-12)*B*C*sin(phi)*sin(theta) ...
% -A*(5.4*sin(5.4*y))*C*sin(5.4*y)*sin(theta)*cos(phi) ; %%f2

%===================== f3 =============================
% f=fc_2d(u)*(-81/16)*(2*(sin(theta)*cos(phi)-0.5)*sin(theta)*(-sin(phi))+2*(sin(theta)*sin(phi)-0.5)*sin(theta)*cos(phi)); %%f3

%===================== f4 =============================
% f = fc_2d(u)*(-81/4)*(2*(sin(theta)*cos(phi)-0.5)*sin(theta)*(-sin(phi))+2*(sin(theta)*sin(phi)-0.5)*sin(theta)*cos(phi)); %%f4

%===================== f5 ============================= 
% f = 0.1*(exp(sin(theta)*cos(phi))*(-sin(phi)*sin(theta)) + exp(sin(theta)*sin(phi)+cos(theta))*cos(phi)*sin(theta)); %%f5

end

function f = der_1y(u)
phi=u(1); theta = u(2);
x = sin(theta)*cos(phi);
y = sin(theta)*sin(phi);
z = cos(theta);

%===================== f1 =============================
A = 0.75*exp(-((9*x-2)^2)/4 - ((9*y-2)^2)/4 - ((9*z-2)^2)/4);
B =  0.75*exp(-((9*x+1)^2)/49 - ((9*y+1)^2)/10 - ((9*z+1)^2)/10);
C = 0.5*exp(-((9*x-7)^2)/4 - ((9*y-3)^2)/4 - ((9*z-5)^2)/4) ;
D = - 0.2*exp (-(9*x-4)^2 - (9*y-7)^2 - (9*z-5)^2);

fx = -0.5*A*(9*x-2) -2/49*B*(9*x+1) -0.5*C*(9*x-7) -2*D*(9*x-4);
fy = -0.5*A*(9*y-2) -0.2*B*(9*y+1) -0.5*C*(9*y-3) -2*D*(9*y-7);
fz = -0.5*A*(9*z-2) -0.2*B*(9*z+1) -0.5*C*(9*z-5) -2*D*(9*z-5);

f = fx*cos(phi)*cos(theta) + fy*cos(theta)*sin(phi) -fz*sin(theta); %%f1
% %===================== f2 =============================
% A = (6+6*(3*x-1)^2)^(-1);
% B = (1.25 + cos(5.4*y));
% C = cos(6*z);
%  
% f=-(6+6*(3*x-1)^2)^(-2)*(36*x-12)*B*C*cos(phi)*cos(theta) ...
% -A*(5.4*sin(5.4*y))*C*sin(5.4*y)*cos(theta)*cos(phi) ...
% +6*A*B*(sin(6*z))*sin(theta); %%f2

%===================== f3 =============================
% f = fc_2d(u)*(-81/16)*(2*(sin(theta)*cos(phi)-0.5)*cos(theta)*cos(phi)+2*(sin(theta)*sin(phi)-0.5)*cos(theta)*sin(phi)+2*(cos(theta)-0.5)*(-sin(theta))); %%f3 

%===================== f4 =============================
% f = fc_2d(u)*(-81/4)*(2*(sin(theta)*cos(phi)-0.5)*cos(theta)*cos(phi)+2*(sin(theta)*sin(phi)-0.5)*cos(theta)*sin(phi)+2*(cos(theta)-0.5)*(-sin(theta))); %%f4 

%===================== f5 =============================
% f = 0.1*(exp(sin(theta)*cos(phi))*cos(theta)*cos(phi) + exp(sin(theta)*sin(phi)+cos(theta))*(cos(theta)*sin(phi)-sin(theta))); %%f5

end 
