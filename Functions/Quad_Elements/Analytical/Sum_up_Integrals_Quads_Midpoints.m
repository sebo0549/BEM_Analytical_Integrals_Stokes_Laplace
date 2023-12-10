function [Ires_ij,Ires_1r]=Sum_up_Integrals_Quads_Midpoints(e1s_t1,e2s_t1,e1s_t2,...
    e2s_t2,e1s_t3,e2s_t3,e1s_t4,e2s_t4,...
    Ixx_t1,Ixy_t1,Iyy_t1,I1r_t1,...
    Ixx_t2,Ixy_t2,Iyy_t2,I1r_t2,...
    Ixx_t3,Ixy_t3,Iyy_t3,I1r_t3,...
    Ixx_t4,Ixy_t4,Iyy_t4,I1r_t4)
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

% perform back transformation for each of the small triangles
Ires_t1=Back_Transformation(e1s_t1,e2s_t1,Ixx_t1,Ixy_t1,Iyy_t1);
Ires_t2=Back_Transformation(e1s_t2,e2s_t2,Ixx_t2,Ixy_t2,Iyy_t2);
Ires_t3=Back_Transformation(e1s_t3,e2s_t3,Ixx_t3,Ixy_t3,Iyy_t3);
Ires_t4=Back_Transformation(e1s_t4,e2s_t4,Ixx_t4,Ixy_t4,Iyy_t4);

% sum up all terms
Ires_ij=Ires_t1+Ires_t2+Ires_t3+Ires_t4;
Ires_1r=I1r_t1+I1r_t2+I1r_t3+I1r_t4;
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