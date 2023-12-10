function [I,F,H]=Reshape_Result_Analytical(Ires_ij,Ires_1r)
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

H(1,1)=Ires_ij(1); 
H(1,2)=Ires_ij(2);
H(1,3)=Ires_ij(3);
H(2,1)=Ires_ij(2);
H(2,2)=Ires_ij(4);
H(2,3)=Ires_ij(5);
H(3,1)=Ires_ij(3);
H(3,2)=Ires_ij(5);
H(3,3)=Ires_ij(6);

F=eye(3)*Ires_1r;
I=H+F;
end