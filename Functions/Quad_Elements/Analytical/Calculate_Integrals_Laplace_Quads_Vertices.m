function [I1r_v1,I1r_v2,I1r_v3,I1r_v4]=Calculate_Integrals_Laplace_Quads_Vertices(P,Con)
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

% vertices of the quads
x1=P(Con(:,1),:);
x2=P(Con(:,2),:);
x3=P(Con(:,3),:);
x4=P(Con(:,4),:);

% shift points
% x1 as source point
x1v1=x2-x1;
x2v1=x3-x1;
x3v1=x4-x1;

% x2 as source point
x1v2=x1-x2;
x2v2=x3-x2;
x3v2=x4-x2;

% x3 as source point
x1v3=x1-x3;
x2v3=x2-x3;
x3v3=x4-x3;

% x4 as source point
x1v4=x1-x4;
x2v4=x2-x4;
x3v4=x3-x4;

% calculate all distances and all angles
[Lq1,phi1v1,phi2v1]=Calculate_Intersection_Points_and_Angles(x3v1,x2v1);
[Lq2,phi3v1,phi4v1]=Calculate_Intersection_Points_and_Angles(x1v1,x2v1);

[Lq3,phi1v2,phi2v2]=Calculate_Intersection_Points_and_Angles(x1v2,x3v2);
[Lq4,phi3v2,phi4v2]=Calculate_Intersection_Points_and_Angles(x2v2,x3v2);

[Lq5,phi1v3,phi2v3]=Calculate_Intersection_Points_and_Angles(x2v3,x1v3);
[Lq6,phi3v3,phi4v3]=Calculate_Intersection_Points_and_Angles(x3v3,x1v3);

[Lq7,phi1v4,phi2v4]=Calculate_Intersection_Points_and_Angles(x3v4,x2v4);
[Lq8,phi3v4,phi4v4]=Calculate_Intersection_Points_and_Angles(x1v4,x2v4);

% calculate the values of the integrals
I1r_v1=Lq1.*log(((1+sin(phi1v1)).*(1+sin(phi2v1)))./(cos(phi1v1).*cos(phi2v1)))+...
       Lq2.*log(((1+sin(phi3v1)).*(1+sin(phi4v1)))./(cos(phi3v1).*cos(phi4v1)));

I1r_v2=Lq3.*log(((1+sin(phi1v2)).*(1+sin(phi2v2)))./(cos(phi1v2).*cos(phi2v2)))+...
       Lq4.*log(((1+sin(phi3v2)).*(1+sin(phi4v2)))./(cos(phi3v2).*cos(phi4v2)));

I1r_v3=Lq5.*log(((1+sin(phi1v3)).*(1+sin(phi2v3)))./(cos(phi1v3).*cos(phi2v3)))+...
       Lq6.*log(((1+sin(phi3v3)).*(1+sin(phi4v3)))./(cos(phi3v3).*cos(phi4v3)));

I1r_v4=Lq7.*log(((1+sin(phi1v4)).*(1+sin(phi2v4)))./(cos(phi1v4).*cos(phi2v4)))+...
       Lq8.*log(((1+sin(phi3v4)).*(1+sin(phi4v4)))./(cos(phi3v4).*cos(phi4v4)));
end