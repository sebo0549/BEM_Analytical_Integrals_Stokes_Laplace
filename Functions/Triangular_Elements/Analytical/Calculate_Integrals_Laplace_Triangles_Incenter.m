function I1r=Calculate_Integrals_Laplace_Triangles_Incenter(P,M_IC,Con)
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

% shift origin off all triangles to the incenters
x1=P(Con(:,1),:)-M_IC;
x2=P(Con(:,2),:)-M_IC;

% calculate all intersection points and all angles
[Lq,phi1,phi2]=Calculate_Intersection_Points_and_Angles(x1,x2);
phi3=pi-phi1-phi2;

% calculate values of the integrals
I1r=2*Lq.*log(((1+sin(phi1)).*(1+sin(phi2)).*(1+sin(phi3)))./(cos(phi1).*cos(phi2).*cos(phi3)));
end