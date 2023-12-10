function [Ires_ij_v1,Ires_1r_v1,...
          Ires_ij_v2,Ires_1r_v2,...
          Ires_ij_v3,Ires_1r_v3,...
          Ires_ij_v4,Ires_1r_v4]=Calculate_Integrals_Stokes_Flow_Quads_Vertices(P,Con)
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
[Lq1v1,phi1v1,phi2v1,q1v1,X1v1]=Calculate_Intersection_Points_and_Angles(x3v1,x2v1);
[Lq2v1,phi3v1,phi4v1,q2v1,X2v1]=Calculate_Intersection_Points_and_Angles(x1v1,x2v1);

[Lq1v2,phi1v2,phi2v2,q1v2,X1v2]=Calculate_Intersection_Points_and_Angles(x1v2,x3v2);
[Lq2v2,phi3v2,phi4v2,q2v2,X2v2]=Calculate_Intersection_Points_and_Angles(x2v2,x3v2);

[Lq1v3,phi1v3,phi2v3,q1v3,X1v3]=Calculate_Intersection_Points_and_Angles(x2v3,x1v3);
[Lq2v3,phi3v3,phi4v3,q2v3,X2v3]=Calculate_Intersection_Points_and_Angles(x3v3,x1v3);

[Lq1v4,phi1v4,phi2v4,q1v4,X1v4]=Calculate_Intersection_Points_and_Angles(x3v4,x2v4);
[Lq2v4,phi3v4,phi4v4,q2v4,X2v4]=Calculate_Intersection_Points_and_Angles(x1v4,x2v4);

% determine matrices for the coordinate transformation
[e1_v1,e2_v1]=Unitary_Transformation(q1v1,X1v1,Lq1v1);
[e3_v1,e4_v1]=Unitary_Transformation(q2v1,X2v1,Lq2v1);

[e1_v2,e2_v2]=Unitary_Transformation(q1v2,X1v2,Lq1v2);
[e3_v2,e4_v2]=Unitary_Transformation(q2v2,X2v2,Lq2v2);

[e1_v3,e2_v3]=Unitary_Transformation(q1v3,X1v3,Lq1v3);
[e3_v3,e4_v3]=Unitary_Transformation(q2v3,X2v3,Lq2v3);

[e1_v4,e2_v4]=Unitary_Transformation(q1v4,X1v4,Lq1v4);
[e3_v4,e4_v4]=Unitary_Transformation(q2v4,X2v4,Lq2v4);

% caclulate integrals for all individual triangles in the new basis
[Ixx1_v1,Ixy1_v1,Iyy1_v1,Irev1_1r_v1]=Integrals_Stokes_Flow(Lq1v1,-phi1v1,phi2v1);
[Ixx2_v1,Ixy2_v1,Iyy2_v1,Irev2_1r_v1]=Integrals_Stokes_Flow(Lq2v1,-phi3v1,phi4v1);

[Ixx1_v2,Ixy1_v2,Iyy1_v2,Irev1_1r_v2]=Integrals_Stokes_Flow(Lq1v2,-phi1v2,phi2v2);
[Ixx2_v2,Ixy2_v2,Iyy2_v2,Irev2_1r_v2]=Integrals_Stokes_Flow(Lq2v2,-phi3v2,phi4v2);

[Ixx1_v3,Ixy1_v3,Iyy1_v3,Irev1_1r_v3]=Integrals_Stokes_Flow(Lq1v3,-phi1v3,phi2v3);
[Ixx2_v3,Ixy2_v3,Iyy2_v3,Irev2_1r_v3]=Integrals_Stokes_Flow(Lq2v3,-phi3v3,phi4v3);

[Ixx1_v4,Ixy1_v4,Iyy1_v4,Irev1_1r_v4]=Integrals_Stokes_Flow(Lq1v4,-phi1v4,phi2v4);
[Ixx2_v4,Ixy2_v4,Iyy2_v4,Irev2_1r_v4]=Integrals_Stokes_Flow(Lq2v4,-phi3v4,phi4v4);

% perform backtransformation and sum up all individual terms
[Ires_ij_v1,Ires_1r_v1,...
    Ires_ij_v2,Ires_1r_v2,...
    Ires_ij_v3,Ires_1r_v3,...
    Ires_ij_v4,Ires_1r_v4]=Sum_up_Integrals_Quads_Vertices(...
                        e1_v1,e2_v1,e3_v1,e4_v1,...
                        e1_v2,e2_v2,e3_v2,e4_v2,...
                        e1_v3,e2_v3,e3_v3,e4_v3,...
                        e1_v4,e2_v4,e3_v4,e4_v4,...
                        Ixx1_v1,Ixy1_v1,Iyy1_v1,Irev1_1r_v1,...
                        Ixx2_v1,Ixy2_v1,Iyy2_v1,Irev2_1r_v1,...
                        Ixx1_v2,Ixy1_v2,Iyy1_v2,Irev1_1r_v2,...
                        Ixx2_v2,Ixy2_v2,Iyy2_v2,Irev2_1r_v2,...
                        Ixx1_v3,Ixy1_v3,Iyy1_v3,Irev1_1r_v3,...
                        Ixx2_v3,Ixy2_v3,Iyy2_v3,Irev2_1r_v3,...
                        Ixx1_v4,Ixy1_v4,Iyy1_v4,Irev1_1r_v4,...
                        Ixx2_v4,Ixy2_v4,Iyy2_v4,Irev2_1r_v4);
end