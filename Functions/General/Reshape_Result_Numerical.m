function [I,F,H]=Reshape_Result_Numerical(Ires_ij,Ires_1r)
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

H=reshape(Ires_ij,3,3);
F=eye(3)*Ires_1r;
I=H+F;
end