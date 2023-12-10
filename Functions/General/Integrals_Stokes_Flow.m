function [Ixx,Ixy,Iyy,I1r]=Integrals_Stokes_Flow(L,uG,oG)
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

% perform integration in the new basis using polar coordinates
Ixx=L.*(sin(oG)-sin(uG));
Ixy=L.*(cos(uG)-cos(oG));
I1r=L.*log( ((1+sin(oG)).*cos(uG))./((1+sin(uG)).*cos(oG)));
Iyy=I1r-Ixx;
end