function [Ires_ij,Ires_1r]=Calculate_Integrals_Stokes_Flow_Triangles_Incenter(P,M_IC,Con)
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

% Shift origin off all triangles to the incenters
x1=P(Con(:,1),:)-M_IC;
x2=P(Con(:,2),:)-M_IC;
x3=P(Con(:,3),:)-M_IC;

% calculate all intersection points and all angles
[Lq,phi1,phi2,q1,X1]=Calculate_Intersection_Points_and_Angles(x1,x2);
[~ ,~   ,phi3,q2,X2]=Calculate_Intersection_Points_and_Angles(x2,x3);
[~ ,~   ,~   ,q3,X3]=Calculate_Intersection_Points_and_Angles(x3,x1);

% determine matrices for the coordinate transformation
[e1s_t1,e2s_t1]=Unitary_Transformation(q1,X1,Lq);
[e1s_t2,e2s_t2]=Unitary_Transformation(q2,X2,Lq);
[e1s_t3,e2s_t3]=Unitary_Transformation(q3,X3,Lq);

% caclulate integrals for all individual triangles in the new basis
[Ixx_t1,Ixy_t1,Iyy_t1,I1r_t1]=Integrals_Stokes_Flow(Lq,-phi1,phi2);
[Ixx_t2,Ixy_t2,Iyy_t2,I1r_t2]=Integrals_Stokes_Flow(Lq,-phi2,phi3);
[Ixx_t3,Ixy_t3,Iyy_t3,I1r_t3]=Integrals_Stokes_Flow(Lq,-phi3,phi1);

% perform backtransformation and sum up all the individual contributions
[Ires_ij,Ires_1r]=Sum_up_Integrals_Triangles_Midpoints(e1s_t1,e2s_t1,e1s_t2,e2s_t2,e1s_t3,e2s_t3,...
    Ixx_t1,Ixy_t1,Iyy_t1,I1r_t1,...
    Ixx_t2,Ixy_t2,Iyy_t2,I1r_t2,...
    Ixx_t3,Ixy_t3,Iyy_t3,I1r_t3);
end