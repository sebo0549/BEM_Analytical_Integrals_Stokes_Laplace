function [Ires_ij,Ires_1r]=Calculate_Integrals_Stokes_Flow_Quads_Midpoints(P,M_Quads,Con)
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

% Shift origin off all quads to the midpoints
x1=P(Con(:,1),:)-M_Quads;
x2=P(Con(:,2),:)-M_Quads;
x3=P(Con(:,3),:)-M_Quads;
x4=P(Con(:,4),:)-M_Quads;

% calculate all distances and all angles
[Lq1,phi1,phi2,q1,X1]=Calculate_Intersection_Points_and_Angles(x1,x2);
[Lq2,phi3,phi4,q2,X2]=Calculate_Intersection_Points_and_Angles(x2,x3);
[Lq3,phi5,phi6,q3,X3]=Calculate_Intersection_Points_and_Angles(x3,x4);
[Lq4,phi7,phi8,q4,X4]=Calculate_Intersection_Points_and_Angles(x4,x1);

% determine matrices for the coordinate transformation
[e1s_t1,e2s_t1]=Unitary_Transformation(q1,X1,Lq1);
[e1s_t2,e2s_t2]=Unitary_Transformation(q2,X2,Lq2);
[e1s_t3,e2s_t3]=Unitary_Transformation(q3,X3,Lq3);
[e1s_t4,e2s_t4]=Unitary_Transformation(q4,X4,Lq4);

% caclulate integrals for all individual triangles in the new basis
[Ixx_t1,Ixy_t1,Iyy_t1,I1r_t1]=Integrals_Stokes_Flow(Lq1,-phi1,phi2);
[Ixx_t2,Ixy_t2,Iyy_t2,I1r_t2]=Integrals_Stokes_Flow(Lq2,-phi3,phi4);
[Ixx_t3,Ixy_t3,Iyy_t3,I1r_t3]=Integrals_Stokes_Flow(Lq3,-phi5,phi6);
[Ixx_t4,Ixy_t4,Iyy_t4,I1r_t4]=Integrals_Stokes_Flow(Lq4,-phi7,phi8);

% perform backtransformation and sum up all the individual contributions
[Ires_ij,Ires_1r]=Sum_up_Integrals_Quads_Midpoints(e1s_t1,e2s_t1,e1s_t2,e2s_t2,e1s_t3,e2s_t3,e1s_t4,e2s_t4,...
    Ixx_t1,Ixy_t1,Iyy_t1,I1r_t1,...
    Ixx_t2,Ixy_t2,Iyy_t2,I1r_t2,...
    Ixx_t3,Ixy_t3,Iyy_t3,I1r_t3,...
    Ixx_t4,Ixy_t4,Iyy_t4,I1r_t4);
end