function [Lq,phi1,phi2,qi,AB]=Calculate_Intersection_Points_and_Angles(x1,x2)
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

% edge vector
AB=x2-x1;

% intersection point
beta=-sum(x1.*AB,2)./sum(AB.*AB,2);
qi=x1+beta.*AB;

% lengths
Lq=sqrt(sum(qi.^2,2));
L1=sqrt(sum(x1.^2,2));  
L2=sqrt(sum(x2.^2,2)); 

% angles
phi1= sign(beta)  .*acos(Lq./L1);
phi2= sign(1-beta).*acos(Lq./L2);
end