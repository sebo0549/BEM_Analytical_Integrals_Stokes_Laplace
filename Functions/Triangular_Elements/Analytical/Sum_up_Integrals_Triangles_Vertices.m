function [Ires_ij_s1,Ires_ij_s2,Ires_ij_s3]=...
    Sum_up_Integrals_Triangles_Vertices(...
    e1s_s1,e2s_s1,e1s_s2,e2s_s2,e1s_s3,e2s_s3,...
    Ixx_s1,Ixy_s1,Iyy_s1,...
    Ixx_s2,Ixy_s2,Iyy_s2,...
    Ixx_s3,Ixy_s3,Iyy_s3)
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

% perform back transformation for each of the small triangles
Ires_ij_s1=Back_Transformation(e1s_s1,e2s_s1,Ixx_s1,Ixy_s1,Iyy_s1);
Ires_ij_s2=Back_Transformation(e1s_s2,e2s_s2,Ixx_s2,Ixy_s2,Iyy_s2);
Ires_ij_s3=Back_Transformation(e1s_s3,e2s_s3,Ixx_s3,Ixy_s3,Iyy_s3);
end

function I_res=Back_Transformation(e1,e2,Ixx_s,Ixy_s,Iyy_s)
I_res=zeros(length(Ixx_s),6);

x_Term1=e1(:,1).*Ixx_s+e2(:,1).*Ixy_s;
x_Term2=e1(:,1).*Ixy_s+e2(:,1).*Iyy_s;
y_Term1=e1(:,2).*Ixx_s+e2(:,2).*Ixy_s;
y_Term2=e1(:,2).*Ixy_s+e2(:,2).*Iyy_s;
z_Term1=e1(:,3).*Ixx_s+e2(:,3).*Ixy_s;
z_Term2=e1(:,3).*Ixy_s+e2(:,3).*Iyy_s;

% Ixx, Ixy, Ixz, Iyy, Iyz, Izz
I_res(:,1)=x_Term1.*e1(:,1) + x_Term2.*e2(:,1);
I_res(:,2)=x_Term1.*e1(:,2) + x_Term2.*e2(:,2);
I_res(:,3)=x_Term1.*e1(:,3) + x_Term2.*e2(:,3);
I_res(:,4)=y_Term1.*e1(:,2) + y_Term2.*e2(:,2);
I_res(:,5)=y_Term1.*e1(:,3) + y_Term2.*e2(:,3);
I_res(:,6)=z_Term1.*e1(:,3) + z_Term2.*e2(:,3);
end