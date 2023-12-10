function [Ires_ij_v1,Ires_1r_v1,...
    Ires_ij_v2,Ires_1r_v2,...
    Ires_ij_v3,Ires_1r_v3]=Calculate_Integrals_Stokes_Flow_Triangles_Vertices(P,Con)
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

% shift origin off all triangles to the incenters
x1=P(Con(:,1),:);
x2=P(Con(:,2),:);
x3=P(Con(:,3),:);

% x1 as source point
x1v1=x2-x1;
x2v1=x3-x1;

% x2 as source point
x1v2=x3-x2;
x2v2=x1-x2;

% x3 as source point
x1v3=x1-x3;
x2v3=x2-x3;

% calculate all intersection points and all angles
[Lqv1,phi1v1,phi2v1,qv1,Xv1]=Calculate_Intersection_Points_and_Angles(x1v1,x2v1);
[Lqv2,phi1v2,phi2v2,qv2,Xv2]=Calculate_Intersection_Points_and_Angles(x1v2,x2v2);
[Lqv3,phi1v3,phi2v3,qv3,Xv3]=Calculate_Intersection_Points_and_Angles(x1v3,x2v3);

% determine matrices for the coordinate transformation
[e1_v1,e2_v1]=Unitary_Transformation(qv1,Xv1,Lqv1);
[e1_v2,e2_v2]=Unitary_Transformation(qv2,Xv2,Lqv2);
[e1_v3,e2_v3]=Unitary_Transformation(qv3,Xv3,Lqv3);

% caclulate integrals for all individual triangles in the new basis
[Ixx_v1,Ixy_v1,Iyy_v1,Ires_1r_v1]=Integrals_Stokes_Flow(Lqv1,-phi1v1,phi2v1);
[Ixx_v2,Ixy_v2,Iyy_v2,Ires_1r_v2]=Integrals_Stokes_Flow(Lqv2,-phi1v2,phi2v2);
[Ixx_v3,Ixy_v3,Iyy_v3,Ires_1r_v3]=Integrals_Stokes_Flow(Lqv3,-phi1v3,phi2v3);

% perform backtransformation
[Ires_ij_v1,Ires_ij_v2,Ires_ij_v3]=...
    Sum_up_Integrals_Triangles_Vertices(...
    e1_v1,e2_v1,e1_v2,e2_v2,e1_v3,e2_v3,...
    Ixx_v1,Ixy_v1,Iyy_v1,...
    Ixx_v2,Ixy_v2,Iyy_v2,...
    Ixx_v3,Ixy_v3,Iyy_v3);
end