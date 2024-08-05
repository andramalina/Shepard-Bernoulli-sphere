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
results_for_rmse= zeros(1, nr_points);
results_for_mae = zeros(1, nr_points);
true_vals = zeros(1, nr_points);
a=shepard_bern_1(nodes(3,:), n, nodes, triangles, miu);
b=fc_3d(nodes(3,:));

%====================================================
for j=1:nr_points
    a=0; b=0;
    a=shepard_bern_1(points(j,:), n, nodes, triangles, miu);
    b=fc_3d(points(j,:));
    results_for_mae(j) = abs(a-b);
    results_for_rmse(j)=a;
    true_vals(j)=b;
end

%====================================================
max_abs_err = max(results_for_mae)
mean_abs_err = mean(results_for_mae)
rmse=(sum((results_for_rmse-true_vals).^2)/nr_points)^1/2

toc
%====================================================

function s = shepard_bern_1(x, n, nodes, triangles, miu)
    s=0;
    numer = sum_ai_all_nodes(x,nodes,miu);
    for i=1:n
        [res_bern,tr] = find_best_triangle_i(x, i, nodes, triangles);
        s = s + (ai(x,i,nodes,miu)/numer)*res_bern;
    end
end

function [result,tr] = find_best_triangle_i(x,i,nodes, triangles)
    all_tr = triangles_with_vertex_i(triangles,i);  
    results = [];
    for j=1:length(all_tr)   
        node1 = nodes(all_tr(j,1),:); 
        node2 = nodes(all_tr(j,2),:);
        node3 = nodes(all_tr(j,3),:);
        
        if det([node1;node2;node3])<0
            det([node1;node2;node3])
            node3=aux;
            node3=node1;
            node1=aux;
            det([node1;node2;node3])
        end
    
        spher_area = spherical_area(node1, node2,node3);
        lambda1 = spherical_area(x,node2,node3)/spher_area;
        lambda2 = spherical_area(node1,x,node3)/spher_area;
        lambda3 = spherical_area(node1,node2,x)/spher_area;
    
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
        results = [results; bern_ord_1(node1, node2, node3, lambda1, lambda2, lambda3)];
    end
    [min_err, index] = min(abs(fc_2d(point)-results));
    result = results(index);
    tr = all_tr(index,:);
end

function tr_vertex_i = triangles_with_vertex_i(triangles,i)
    [row,~]=find(triangles==i);
    tr_vertex_i = triangles(row,:);
end

function f =  bern_ord_1(node1, node2, node3, lambda1, lambda2, lambda3)
    f = 0;
    f = lambda1*fc_2d(node1) +  lambda2*fc_2d(node2)+ lambda3*fc_2d(node3);
end

function r = ai(x,i,nodes, miu)
    r = (1/geodesic_dist(x,nodes(i,:)))^miu;
end

function s = sum_ai_all_nodes(x,nodes,miu)
    s=0;
    n=length(nodes);
    for i = 1:n
        s= s+ geodesic_dist(x,nodes(i,:))^(-miu);
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