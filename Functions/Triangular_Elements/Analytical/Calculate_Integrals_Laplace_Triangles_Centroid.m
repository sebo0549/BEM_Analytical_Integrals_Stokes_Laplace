function I1r=Calculate_Integrals_Laplace_Triangles_Centroid(P,M_Centroid,Con)
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

% shift origin off all triangles to the incenters
x1=P(Con(:,1),:)-M_Centroid;
x2=P(Con(:,2),:)-M_Centroid;
x3=P(Con(:,3),:)-M_Centroid;

% calculate all intersection points and all angles
[Lq1,phi1,phi2]=Calculate_Intersection_Points_and_Angles(x1,x2);
[Lq2,phi3,phi4]=Calculate_Intersection_Points_and_Angles(x2,x3);
[Lq3,phi5,phi6]=Calculate_Intersection_Points_and_Angles(x3,x1);

% calculate values of the integral
I1r=Lq1.*log(((1+sin(phi1)).*(1+sin(phi2)))./(cos(phi1).*cos(phi2)))+...
    Lq2.*log(((1+sin(phi3)).*(1+sin(phi4)))./(cos(phi3).*cos(phi4)))+...
    Lq3.*log(((1+sin(phi5)).*(1+sin(phi6)))./(cos(phi5).*cos(phi6)));
end