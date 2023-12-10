function [Ires_ij_s1,Ires_1r_s1,...
    Ires_ij_s2,Ires_1r_s2,...
    Ires_ij_s3,Ires_1r_s3]=Calculate_Integrals_Numerical_Triangles_Vertices(P,abs_tol)
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

warning('off','MATLAB:quad2d:minRectSizePass')
warning('off','MATLAB:quad2d:maxFunEvalsFail')
warning('off','MATLAB:quad2d:maxFunEvalsPass')

% Shift origin  to the midpoint
Ps1=P-P(1,:);
Ps2=P-P(2,:);
Ps3=P-P(3,:);

% s1:
As1_t=0.5*norm(cross(Ps1(2,:),Ps1(3,:)));
rs1_t_x=@(eta,psi) Ps1(2,1)*eta+Ps1(3,1)*psi;
rs1_t_y=@(eta,psi) Ps1(2,2)*eta+Ps1(3,2)*psi;
rs1_t_z=@(eta,psi) Ps1(2,3)*eta+Ps1(3,3)*psi;

% s2:
As2_t=0.5*norm(cross(Ps2(1,:),Ps2(3,:)));
rs2_t_x=@(eta,psi) Ps2(1,1)*eta+Ps2(3,1)*psi;
rs2_t_y=@(eta,psi) Ps2(1,2)*eta+Ps2(3,2)*psi;
rs2_t_z=@(eta,psi) Ps2(1,3)*eta+Ps2(3,3)*psi;

% s3:
As3_t=0.5*norm(cross(Ps3(1,:),Ps3(2,:)));
rs3_t_x=@(eta,psi) Ps3(1,1)*eta+Ps3(2,1)*psi;
rs3_t_y=@(eta,psi) Ps3(1,2)*eta+Ps3(2,2)*psi;
rs3_t_z=@(eta,psi) Ps3(1,3)*eta+Ps3(2,3)*psi;

% 1/r
f_1r_s1=@(eta,psi) 2*As1_t*1./(rs1_t_x(eta,psi).^2+rs1_t_y(eta,psi).^2+rs1_t_z(eta,psi).^2).^(1/2);
f_1r_s2=@(eta,psi) 2*As2_t*1./(rs2_t_x(eta,psi).^2+rs2_t_y(eta,psi).^2+rs2_t_z(eta,psi).^2).^(1/2);
f_1r_s3=@(eta,psi) 2*As3_t*1./(rs3_t_x(eta,psi).^2+rs3_t_y(eta,psi).^2+rs3_t_z(eta,psi).^2).^(1/2);

% xx
f_xx_s1=@(eta,psi) 2*As1_t*(rs1_t_x(eta,psi).*rs1_t_x(eta,psi))./(rs1_t_x(eta,psi).^2+rs1_t_y(eta,psi).^2+rs1_t_z(eta,psi).^2).^(3/2);
f_xx_s2=@(eta,psi) 2*As2_t*(rs2_t_x(eta,psi).*rs2_t_x(eta,psi))./(rs2_t_x(eta,psi).^2+rs2_t_y(eta,psi).^2+rs2_t_z(eta,psi).^2).^(3/2);
f_xx_s3=@(eta,psi) 2*As3_t*(rs3_t_x(eta,psi).*rs3_t_x(eta,psi))./(rs3_t_x(eta,psi).^2+rs3_t_y(eta,psi).^2+rs3_t_z(eta,psi).^2).^(3/2);

% xy
f_xy_s1=@(eta,psi) 2*As1_t*(rs1_t_x(eta,psi).*rs1_t_y(eta,psi))./(rs1_t_x(eta,psi).^2+rs1_t_y(eta,psi).^2+rs1_t_z(eta,psi).^2).^(3/2);
f_xy_s2=@(eta,psi) 2*As2_t*(rs2_t_x(eta,psi).*rs2_t_y(eta,psi))./(rs2_t_x(eta,psi).^2+rs2_t_y(eta,psi).^2+rs2_t_z(eta,psi).^2).^(3/2);
f_xy_s3=@(eta,psi) 2*As3_t*(rs3_t_x(eta,psi).*rs3_t_y(eta,psi))./(rs3_t_x(eta,psi).^2+rs3_t_y(eta,psi).^2+rs3_t_z(eta,psi).^2).^(3/2);

% xz
f_xz_s1=@(eta,psi) 2*As1_t*(rs1_t_x(eta,psi).*rs1_t_z(eta,psi))./(rs1_t_x(eta,psi).^2+rs1_t_y(eta,psi).^2+rs1_t_z(eta,psi).^2).^(3/2);
f_xz_s2=@(eta,psi) 2*As2_t*(rs2_t_x(eta,psi).*rs2_t_z(eta,psi))./(rs2_t_x(eta,psi).^2+rs2_t_y(eta,psi).^2+rs2_t_z(eta,psi).^2).^(3/2);
f_xz_s3=@(eta,psi) 2*As3_t*(rs3_t_x(eta,psi).*rs3_t_z(eta,psi))./(rs3_t_x(eta,psi).^2+rs3_t_y(eta,psi).^2+rs3_t_z(eta,psi).^2).^(3/2);

% yy
f_yy_s1=@(eta,psi) 2*As1_t*(rs1_t_y(eta,psi).*rs1_t_y(eta,psi))./(rs1_t_x(eta,psi).^2+rs1_t_y(eta,psi).^2+rs1_t_z(eta,psi).^2).^(3/2);
f_yy_s2=@(eta,psi) 2*As2_t*(rs2_t_y(eta,psi).*rs2_t_y(eta,psi))./(rs2_t_x(eta,psi).^2+rs2_t_y(eta,psi).^2+rs2_t_z(eta,psi).^2).^(3/2);
f_yy_s3=@(eta,psi) 2*As3_t*(rs3_t_y(eta,psi).*rs3_t_y(eta,psi))./(rs3_t_x(eta,psi).^2+rs3_t_y(eta,psi).^2+rs3_t_z(eta,psi).^2).^(3/2);

% zy
f_zy_s1=@(eta,psi) 2*As1_t*(rs1_t_z(eta,psi).*rs1_t_y(eta,psi))./(rs1_t_x(eta,psi).^2+rs1_t_y(eta,psi).^2+rs1_t_z(eta,psi).^2).^(3/2);
f_zy_s2=@(eta,psi) 2*As2_t*(rs2_t_z(eta,psi).*rs2_t_y(eta,psi))./(rs2_t_x(eta,psi).^2+rs2_t_y(eta,psi).^2+rs2_t_z(eta,psi).^2).^(3/2);
f_zy_s3=@(eta,psi) 2*As3_t*(rs3_t_z(eta,psi).*rs3_t_y(eta,psi))./(rs3_t_x(eta,psi).^2+rs3_t_y(eta,psi).^2+rs3_t_z(eta,psi).^2).^(3/2);

% zz
f_zz_s1=@(eta,psi) 2*As1_t*(rs1_t_z(eta,psi).*rs1_t_z(eta,psi))./(rs1_t_x(eta,psi).^2+rs1_t_y(eta,psi).^2+rs1_t_z(eta,psi).^2).^(3/2);
f_zz_s2=@(eta,psi) 2*As2_t*(rs2_t_z(eta,psi).*rs2_t_z(eta,psi))./(rs2_t_x(eta,psi).^2+rs2_t_y(eta,psi).^2+rs2_t_z(eta,psi).^2).^(3/2);
f_zz_s3=@(eta,psi) 2*As3_t*(rs3_t_z(eta,psi).*rs3_t_z(eta,psi))./(rs3_t_x(eta,psi).^2+rs3_t_y(eta,psi).^2+rs3_t_z(eta,psi).^2).^(3/2);

% perform numerical integration
ymax = @(x) 1 - x;

Ires_1r_s1=quad2d(f_1r_s1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_1r_s2=quad2d(f_1r_s2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_1r_s3=quad2d(f_1r_s3,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);

Ires_ij_s1(1)=quad2d(f_xx_s1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_s2(1)=quad2d(f_xx_s2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_s3(1)=quad2d(f_xx_s3,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);

Ires_ij_s1(2)=quad2d(f_xy_s1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_s2(2)=quad2d(f_xy_s2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_s3(2)=quad2d(f_xy_s3,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);

Ires_ij_s1(3)=quad2d(f_xz_s1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_s2(3)=quad2d(f_xz_s2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_s3(3)=quad2d(f_xz_s3,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);

Ires_ij_s1(5)=quad2d(f_yy_s1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_s2(5)=quad2d(f_yy_s2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_s3(5)=quad2d(f_yy_s3,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);

Ires_ij_s1(8)=quad2d(f_zy_s1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_s2(8)=quad2d(f_zy_s2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_s3(8)=quad2d(f_zy_s3,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);

Ires_ij_s1(9)=quad2d(f_zz_s1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_s2(9)=quad2d(f_zz_s2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_s3(9)=quad2d(f_zz_s3,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);

Ires_ij_s1(4)=Ires_ij_s1(2);
Ires_ij_s1(6)=Ires_ij_s1(8);
Ires_ij_s1(7)=Ires_ij_s1(3);

Ires_ij_s2(4)=Ires_ij_s2(2);
Ires_ij_s2(6)=Ires_ij_s2(8);
Ires_ij_s2(7)=Ires_ij_s2(3);

Ires_ij_s3(4)=Ires_ij_s3(2);
Ires_ij_s3(6)=Ires_ij_s3(8);
Ires_ij_s3(7)=Ires_ij_s3(3);
end