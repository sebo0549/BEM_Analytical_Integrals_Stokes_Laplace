function [Ires_ij_s1,Ires_1r_s1,...
    Ires_ij_s2,Ires_1r_s2,...
    Ires_ij_s3,Ires_1r_s3,...
    Ires_ij_s4,Ires_1r_s4]=...
    Sum_up_Integrals_Quads_Vertices(...
    e1_s1,e2_s1,e3_s1,e4_s1,...
    e1_s2,e2_s2,e3_s2,e4_s2,...
    e1_s3,e2_s3,e3_s3,e4_s3,...
    e1_s4,e2_s4,e3_s4,e4_s4,...
    Ixx1_s1,Ixy1_s1,Iyy1_s1,Ires1_1r_s1,...
    Ixx2_s1,Ixy2_s1,Iyy2_s1,Ires2_1r_s1,...
    Ixx1_s2,Ixy1_s2,Iyy1_s2,Ires1_1r_s2,...
    Ixx2_s2,Ixy2_s2,Iyy2_s2,Ires2_1r_s2,...
    Ixx1_s3,Ixy1_s3,Iyy1_s3,Ires1_1r_s3,...
    Ixx2_s3,Ixy2_s3,Iyy2_s3,Ires2_1r_s3,...
    Ixx1_s4,Ixy1_s4,Iyy1_s4,Ires1_1r_s4,...
    Ixx2_s4,Ixy2_s4,Iyy2_s4,Ires2_1r_s4)
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

% perform back transformation for each of the small triangles
Ires1_s1=Back_Transformation(e1_s1,e2_s1,Ixx1_s1,Ixy1_s1,Iyy1_s1);
Ires2_s1=Back_Transformation(e3_s1,e4_s1,Ixx2_s1,Ixy2_s1,Iyy2_s1);

Ires1_s2=Back_Transformation(e1_s2,e2_s2,Ixx1_s2,Ixy1_s2,Iyy1_s2);
Ires2_s2=Back_Transformation(e3_s2,e4_s2,Ixx2_s2,Ixy2_s2,Iyy2_s2);

Ires1_s3=Back_Transformation(e1_s3,e2_s3,Ixx1_s3,Ixy1_s3,Iyy1_s3);
Ires2_s3=Back_Transformation(e3_s3,e4_s3,Ixx2_s3,Ixy2_s3,Iyy2_s3);

Ires1_s4=Back_Transformation(e1_s4,e2_s4,Ixx1_s4,Ixy1_s4,Iyy1_s4);
Ires2_s4=Back_Transformation(e3_s4,e4_s4,Ixx2_s4,Ixy2_s4,Iyy2_s4);

% sum up all terms
Ires_ij_s1=Ires1_s1+Ires2_s1;
Ires_1r_s1=Ires1_1r_s1+Ires2_1r_s1;

Ires_ij_s2=Ires1_s2+Ires2_s2;
Ires_1r_s2=Ires1_1r_s2+Ires2_1r_s2;

Ires_ij_s3=Ires1_s3+Ires2_s3;
Ires_1r_s3=Ires1_1r_s3+Ires2_1r_s3;

Ires_ij_s4=Ires1_s4+Ires2_s4;
Ires_1r_s4=Ires1_1r_s4+Ires2_1r_s4;
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