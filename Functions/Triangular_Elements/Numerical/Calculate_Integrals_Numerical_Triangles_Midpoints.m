function [Ires_ij,Ires_1r]=Calculate_Integrals_Numerical_Triangles_Midpoints(P,MP,abs_tol)
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

warning('off','MATLAB:quad2d:minRectSizePass')
warning('off','MATLAB:quad2d:maxFunEvalsFail')
warning('off','MATLAB:quad2d:maxFunEvalsPass')

% Shift origin  to the midpoint
PM=P-MP;

% T1:
A_t1=0.5*norm(cross(PM(1,:),PM(2,:)));
r_t1_x=@(eta,psi) PM(1,1)*eta+PM(2,1)*psi;
r_t1_y=@(eta,psi) PM(1,2)*eta+PM(2,2)*psi;
r_t1_z=@(eta,psi) PM(1,3)*eta+PM(2,3)*psi;

% T2:
A_t2=0.5*norm(cross(PM(2,:),PM(3,:)));
r_t2_x=@(eta,psi) PM(2,1)*eta+PM(3,1)*psi;
r_t2_y=@(eta,psi) PM(2,2)*eta+PM(3,2)*psi;
r_t2_z=@(eta,psi) PM(2,3)*eta+PM(3,3)*psi;

% T3:
A_t3=0.5*norm(cross(PM(3,:),PM(1,:)));
r_t3_x=@(eta,psi) PM(1,1)*eta+PM(3,1)*psi;
r_t3_y=@(eta,psi) PM(1,2)*eta+PM(3,2)*psi;
r_t3_z=@(eta,psi) PM(1,3)*eta+PM(3,3)*psi;

% 1/r
f_1r_t1=@(eta,psi) 2*A_t1*1./(r_t1_x(eta,psi).^2+r_t1_y(eta,psi).^2+r_t1_z(eta,psi).^2).^(1/2);
f_1r_t2=@(eta,psi) 2*A_t2*1./(r_t2_x(eta,psi).^2+r_t2_y(eta,psi).^2+r_t2_z(eta,psi).^2).^(1/2);
f_1r_t3=@(eta,psi) 2*A_t3*1./(r_t3_x(eta,psi).^2+r_t3_y(eta,psi).^2+r_t3_z(eta,psi).^2).^(1/2);

% xx
f_xx_t1=@(eta,psi) 2*A_t1*(r_t1_x(eta,psi).*r_t1_x(eta,psi))./(r_t1_x(eta,psi).^2+r_t1_y(eta,psi).^2+r_t1_z(eta,psi).^2).^(3/2);
f_xx_t2=@(eta,psi) 2*A_t2*(r_t2_x(eta,psi).*r_t2_x(eta,psi))./(r_t2_x(eta,psi).^2+r_t2_y(eta,psi).^2+r_t2_z(eta,psi).^2).^(3/2);
f_xx_t3=@(eta,psi) 2*A_t3*(r_t3_x(eta,psi).*r_t3_x(eta,psi))./(r_t3_x(eta,psi).^2+r_t3_y(eta,psi).^2+r_t3_z(eta,psi).^2).^(3/2);

% xy
f_xy_t1=@(eta,psi) 2*A_t1*(r_t1_x(eta,psi).*r_t1_y(eta,psi))./(r_t1_x(eta,psi).^2+r_t1_y(eta,psi).^2+r_t1_z(eta,psi).^2).^(3/2);
f_xy_t2=@(eta,psi) 2*A_t2*(r_t2_x(eta,psi).*r_t2_y(eta,psi))./(r_t2_x(eta,psi).^2+r_t2_y(eta,psi).^2+r_t2_z(eta,psi).^2).^(3/2);
f_xy_t3=@(eta,psi) 2*A_t3*(r_t3_x(eta,psi).*r_t3_y(eta,psi))./(r_t3_x(eta,psi).^2+r_t3_y(eta,psi).^2+r_t3_z(eta,psi).^2).^(3/2);

% xz
f_xz_t1=@(eta,psi) 2*A_t1*(r_t1_x(eta,psi).*r_t1_z(eta,psi))./(r_t1_x(eta,psi).^2+r_t1_y(eta,psi).^2+r_t1_z(eta,psi).^2).^(3/2);
f_xz_t2=@(eta,psi) 2*A_t2*(r_t2_x(eta,psi).*r_t2_z(eta,psi))./(r_t2_x(eta,psi).^2+r_t2_y(eta,psi).^2+r_t2_z(eta,psi).^2).^(3/2);
f_xz_t3=@(eta,psi) 2*A_t3*(r_t3_x(eta,psi).*r_t3_z(eta,psi))./(r_t3_x(eta,psi).^2+r_t3_y(eta,psi).^2+r_t3_z(eta,psi).^2).^(3/2);

% yy
f_yy_t1=@(eta,psi) 2*A_t1*(r_t1_y(eta,psi).*r_t1_y(eta,psi))./(r_t1_x(eta,psi).^2+r_t1_y(eta,psi).^2+r_t1_z(eta,psi).^2).^(3/2);
f_yy_t2=@(eta,psi) 2*A_t2*(r_t2_y(eta,psi).*r_t2_y(eta,psi))./(r_t2_x(eta,psi).^2+r_t2_y(eta,psi).^2+r_t2_z(eta,psi).^2).^(3/2);
f_yy_t3=@(eta,psi) 2*A_t3*(r_t3_y(eta,psi).*r_t3_y(eta,psi))./(r_t3_x(eta,psi).^2+r_t3_y(eta,psi).^2+r_t3_z(eta,psi).^2).^(3/2);

% zy
f_zy_t1=@(eta,psi) 2*A_t1*(r_t1_z(eta,psi).*r_t1_y(eta,psi))./(r_t1_x(eta,psi).^2+r_t1_y(eta,psi).^2+r_t1_z(eta,psi).^2).^(3/2);
f_zy_t2=@(eta,psi) 2*A_t2*(r_t2_z(eta,psi).*r_t2_y(eta,psi))./(r_t2_x(eta,psi).^2+r_t2_y(eta,psi).^2+r_t2_z(eta,psi).^2).^(3/2);
f_zy_t3=@(eta,psi) 2*A_t3*(r_t3_z(eta,psi).*r_t3_y(eta,psi))./(r_t3_x(eta,psi).^2+r_t3_y(eta,psi).^2+r_t3_z(eta,psi).^2).^(3/2);

% zz
f_zz_t1=@(eta,psi) 2*A_t1*(r_t1_z(eta,psi).*r_t1_z(eta,psi))./(r_t1_x(eta,psi).^2+r_t1_y(eta,psi).^2+r_t1_z(eta,psi).^2).^(3/2);
f_zz_t2=@(eta,psi) 2*A_t2*(r_t2_z(eta,psi).*r_t2_z(eta,psi))./(r_t2_x(eta,psi).^2+r_t2_y(eta,psi).^2+r_t2_z(eta,psi).^2).^(3/2);
f_zz_t3=@(eta,psi) 2*A_t3*(r_t3_z(eta,psi).*r_t3_z(eta,psi))./(r_t3_x(eta,psi).^2+r_t3_y(eta,psi).^2+r_t3_z(eta,psi).^2).^(3/2);

% perform numerical integration
ymax = @(x) 1 - x;

I1r_t1=quad2d(f_1r_t1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
I1r_t2=quad2d(f_1r_t2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
I1r_t3=quad2d(f_1r_t3,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_1r =I1r_t1+I1r_t2+I1r_t3;

Ixx_t1=quad2d(f_xx_t1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ixx_t2=quad2d(f_xx_t2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ixx_t3=quad2d(f_xx_t3,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij(1) =Ixx_t1+Ixx_t2+Ixx_t3;

Ixy_t1=quad2d(f_xy_t1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ixy_t2=quad2d(f_xy_t2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ixy_t3=quad2d(f_xy_t3,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij(2) =Ixy_t1+Ixy_t2+Ixy_t3;

Ixz_t1=quad2d(f_xz_t1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ixz_t2=quad2d(f_xz_t2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ixz_t3=quad2d(f_xz_t3,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij(3) =Ixz_t1+Ixz_t2+Ixz_t3;

Iyy_t1=quad2d(f_yy_t1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Iyy_t2=quad2d(f_yy_t2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Iyy_t3=quad2d(f_yy_t3,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij(5) =Iyy_t1+Iyy_t2+Iyy_t3;

Izy_t1=quad2d(f_zy_t1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Izy_t2=quad2d(f_zy_t2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Izy_t3=quad2d(f_zy_t3,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij(8) =Izy_t1+Izy_t2+Izy_t3;

Izz_t1=quad2d(f_zz_t1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Izz_t2=quad2d(f_zz_t2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Izz_t3=quad2d(f_zz_t3,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij(9) =Izz_t1+Izz_t2+Izz_t3;

Ires_ij(4)=Ires_ij(2);
Ires_ij(6)=Ires_ij(8);
Ires_ij(7)=Ires_ij(3);
end