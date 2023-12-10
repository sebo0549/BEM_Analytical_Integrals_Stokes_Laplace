function I1r=Calculate_Integrals_Laplace_Quads_Midpoints(P,M_Quads,Con)
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

% Shift origin off all quads to the midpoints
x1=P(Con(:,1),:)-M_Quads;
x2=P(Con(:,2),:)-M_Quads;
x3=P(Con(:,3),:)-M_Quads;
x4=P(Con(:,4),:)-M_Quads;

% calculate all intersection points and all angles
[Lq1,phi1,phi2]=Calculate_Intersection_Points_and_Angles(x1,x2);
[Lq2,phi3,phi4]=Calculate_Intersection_Points_and_Angles(x2,x3);
[Lq3,phi5,phi6]=Calculate_Intersection_Points_and_Angles(x3,x4);
[Lq4,phi7,phi8]=Calculate_Intersection_Points_and_Angles(x4,x1);

% calculate values of the integrals
I1r=Lq1.*log(((1+sin(phi1)).*(1+sin(phi2)))./(cos(phi1).*cos(phi2)))+...
    Lq2.*log(((1+sin(phi3)).*(1+sin(phi4)))./(cos(phi3).*cos(phi4)))+...
    Lq3.*log(((1+sin(phi5)).*(1+sin(phi6)))./(cos(phi5).*cos(phi6)))+...
    Lq4.*log(((1+sin(phi7)).*(1+sin(phi8)))./(cos(phi7).*cos(phi8)));
end