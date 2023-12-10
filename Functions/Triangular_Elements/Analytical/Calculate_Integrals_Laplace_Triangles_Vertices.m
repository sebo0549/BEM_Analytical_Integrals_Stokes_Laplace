function [I1r_v1,I1r_v2,I1r_v3]=Calculate_Integrals_Laplace_Triangles_Vertices(P,Con)
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

% vertices of the triangles
x1=P(Con(:,1),:);
x2=P(Con(:,2),:);
x3=P(Con(:,3),:);

% shift points
% x1 as source point
x1v1=x2-x1;
x1v2=x3-x1;

% x2 as source point
x2v1=x3-x2;
x2v2=x1-x2;

% x3 as source point
x3v1=x1-x3;
x3v2=x2-x3;

% calculate all distances and all angles
[Lq1,phi1v1,phi2v1]=Calculate_Intersection_Points_and_Angles(x1v1,x1v2);
[Lq2,phi1v2,phi2v2]=Calculate_Intersection_Points_and_Angles(x2v1,x2v2);
[Lq3,phi1v3,phi2v3]=Calculate_Intersection_Points_and_Angles(x3v1,x3v2);

% calculate the values of the integrals
I1r_v1=Lq1.*log(((1+sin(phi1v1)).*(1+sin(phi2v1)))./(cos(phi1v1).*cos(phi2v1)));
I1r_v2=Lq2.*log(((1+sin(phi1v2)).*(1+sin(phi2v2)))./(cos(phi1v2).*cos(phi2v2)));
I1r_v3=Lq3.*log(((1+sin(phi1v3)).*(1+sin(phi2v3)))./(cos(phi1v3).*cos(phi2v3)));
end