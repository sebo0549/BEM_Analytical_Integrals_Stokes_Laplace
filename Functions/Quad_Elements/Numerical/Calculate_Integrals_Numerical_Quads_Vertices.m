function [Ires_ij_v1,Ires_1r_v1,...
    Ires_ij_v2,Ires_1r_v2,...
    Ires_ij_v3,Ires_1r_v3,...
    Ires_ij_v4,Ires_1r_v4]=Calculate_Integrals_Numerical_Quads_Vertices(P,abs_tol)
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

warning('off','MATLAB:quad2d:minRectSizePass')
warning('off','MATLAB:quad2d:maxFunEvalsFail')
warning('off','MATLAB:quad2d:maxFunEvalsPass')

% Shift origin  to the midpoint
Pv1=P-P(1,:);
Pv2=P-P(2,:);
Pv3=P-P(3,:);
Pv4=P-P(4,:);

% v1:
Av1_t1=0.5*norm(cross(Pv1(2,:),Pv1(3,:)));
rv1_t1_x=@(eta,psi) Pv1(2,1)*eta+Pv1(3,1)*psi;
rv1_t1_y=@(eta,psi) Pv1(2,2)*eta+Pv1(3,2)*psi;
rv1_t1_z=@(eta,psi) Pv1(2,3)*eta+Pv1(3,3)*psi;

Av1_t2=0.5*norm(cross(Pv1(4,:),Pv1(3,:)));
rv1_t2_x=@(eta,psi) Pv1(4,1)*eta+Pv1(3,1)*psi;
rv1_t2_y=@(eta,psi) Pv1(4,2)*eta+Pv1(3,2)*psi;
rv1_t2_z=@(eta,psi) Pv1(4,3)*eta+Pv1(3,3)*psi;


% v2:
Av2_t1=0.5*norm(cross(Pv2(1,:),Pv2(4,:)));
rv2_t1_x=@(eta,psi) Pv2(1,1)*eta+Pv2(4,1)*psi;
rv2_t1_y=@(eta,psi) Pv2(1,2)*eta+Pv2(4,2)*psi;
rv2_t1_z=@(eta,psi) Pv2(1,3)*eta+Pv2(4,3)*psi;

Av2_t2=0.5*norm(cross(Pv2(3,:),Pv2(4,:)));
rv2_t2_x=@(eta,psi) Pv2(3,1)*eta+Pv2(4,1)*psi;
rv2_t2_y=@(eta,psi) Pv2(3,2)*eta+Pv2(4,2)*psi;
rv2_t2_z=@(eta,psi) Pv2(3,3)*eta+Pv2(4,3)*psi;

% v3:
Av3_t1=0.5*norm(cross(Pv3(2,:),Pv3(1,:)));
rv3_t1_x=@(eta,psi) Pv3(2,1)*eta+Pv3(1,1)*psi;
rv3_t1_y=@(eta,psi) Pv3(2,2)*eta+Pv3(1,2)*psi;
rv3_t1_z=@(eta,psi) Pv3(2,3)*eta+Pv3(1,3)*psi;

Av3_t2=0.5*norm(cross(Pv3(4,:),Pv3(1,:)));
rv3_t2_x=@(eta,psi) Pv3(4,1)*eta+Pv3(1,1)*psi;
rv3_t2_y=@(eta,psi) Pv3(4,2)*eta+Pv3(1,2)*psi;
rv3_t2_z=@(eta,psi) Pv3(4,3)*eta+Pv3(1,3)*psi;

% v4:
Av4_t1=0.5*norm(cross(Pv4(1,:),Pv4(2,:)));
rv4_t1_x=@(eta,psi) Pv4(1,1)*eta+Pv4(2,1)*psi;
rv4_t1_y=@(eta,psi) Pv4(1,2)*eta+Pv4(2,2)*psi;
rv4_t1_z=@(eta,psi) Pv4(1,3)*eta+Pv4(2,3)*psi;

Av4_t2=0.5*norm(cross(Pv4(3,:),Pv4(2,:)));
rv4_t2_x=@(eta,psi) Pv4(3,1)*eta+Pv4(2,1)*psi;
rv4_t2_y=@(eta,psi) Pv4(3,2)*eta+Pv4(2,2)*psi;
rv4_t2_z=@(eta,psi) Pv4(3,3)*eta+Pv4(2,3)*psi;


% 1/r
f_1r_v1_1=@(eta,psi) 2*Av1_t1*1./(rv1_t1_x(eta,psi).^2+rv1_t1_y(eta,psi).^2+rv1_t1_z(eta,psi).^2).^(1/2);
f_1r_v1_2=@(eta,psi) 2*Av1_t2*1./(rv1_t2_x(eta,psi).^2+rv1_t2_y(eta,psi).^2+rv1_t2_z(eta,psi).^2).^(1/2);
f_1r_v2_1=@(eta,psi) 2*Av2_t1*1./(rv2_t1_x(eta,psi).^2+rv2_t1_y(eta,psi).^2+rv2_t1_z(eta,psi).^2).^(1/2);
f_1r_v2_2=@(eta,psi) 2*Av2_t2*1./(rv2_t2_x(eta,psi).^2+rv2_t2_y(eta,psi).^2+rv2_t2_z(eta,psi).^2).^(1/2);
f_1r_v3_1=@(eta,psi) 2*Av3_t1*1./(rv3_t1_x(eta,psi).^2+rv3_t1_y(eta,psi).^2+rv3_t1_z(eta,psi).^2).^(1/2);
f_1r_v3_2=@(eta,psi) 2*Av3_t2*1./(rv3_t2_x(eta,psi).^2+rv3_t2_y(eta,psi).^2+rv3_t2_z(eta,psi).^2).^(1/2);
f_1r_v4_1=@(eta,psi) 2*Av4_t1*1./(rv4_t1_x(eta,psi).^2+rv4_t1_y(eta,psi).^2+rv4_t1_z(eta,psi).^2).^(1/2);
f_1r_v4_2=@(eta,psi) 2*Av4_t2*1./(rv4_t2_x(eta,psi).^2+rv4_t2_y(eta,psi).^2+rv4_t2_z(eta,psi).^2).^(1/2);

% xx
f_xx_v1_1=@(eta,psi) 2*Av1_t1*(rv1_t1_x(eta,psi).*rv1_t1_x(eta,psi))./(rv1_t1_x(eta,psi).^2+rv1_t1_y(eta,psi).^2+rv1_t1_z(eta,psi).^2).^(3/2);
f_xx_v1_2=@(eta,psi) 2*Av1_t2*(rv1_t2_x(eta,psi).*rv1_t2_x(eta,psi))./(rv1_t2_x(eta,psi).^2+rv1_t2_y(eta,psi).^2+rv1_t2_z(eta,psi).^2).^(3/2);
f_xx_v2_1=@(eta,psi) 2*Av2_t1*(rv2_t1_x(eta,psi).*rv2_t1_x(eta,psi))./(rv2_t1_x(eta,psi).^2+rv2_t1_y(eta,psi).^2+rv2_t1_z(eta,psi).^2).^(3/2);
f_xx_v2_2=@(eta,psi) 2*Av2_t2*(rv2_t2_x(eta,psi).*rv2_t2_x(eta,psi))./(rv2_t2_x(eta,psi).^2+rv2_t2_y(eta,psi).^2+rv2_t2_z(eta,psi).^2).^(3/2);
f_xx_v3_1=@(eta,psi) 2*Av3_t1*(rv3_t1_x(eta,psi).*rv3_t1_x(eta,psi))./(rv3_t1_x(eta,psi).^2+rv3_t1_y(eta,psi).^2+rv3_t1_z(eta,psi).^2).^(3/2);
f_xx_v3_2=@(eta,psi) 2*Av3_t2*(rv3_t2_x(eta,psi).*rv3_t2_x(eta,psi))./(rv3_t2_x(eta,psi).^2+rv3_t2_y(eta,psi).^2+rv3_t2_z(eta,psi).^2).^(3/2);
f_xx_v4_1=@(eta,psi) 2*Av4_t1*(rv4_t1_x(eta,psi).*rv4_t1_x(eta,psi))./(rv4_t1_x(eta,psi).^2+rv4_t1_y(eta,psi).^2+rv4_t1_z(eta,psi).^2).^(3/2);
f_xx_v4_2=@(eta,psi) 2*Av4_t2*(rv4_t2_x(eta,psi).*rv4_t2_x(eta,psi))./(rv4_t2_x(eta,psi).^2+rv4_t2_y(eta,psi).^2+rv4_t2_z(eta,psi).^2).^(3/2);

% xy
f_xy_v1_1=@(eta,psi) 2*Av1_t1*(rv1_t1_x(eta,psi).*rv1_t1_y(eta,psi))./(rv1_t1_x(eta,psi).^2+rv1_t1_y(eta,psi).^2+rv1_t1_z(eta,psi).^2).^(3/2);
f_xy_v1_2=@(eta,psi) 2*Av1_t2*(rv1_t2_x(eta,psi).*rv1_t2_y(eta,psi))./(rv1_t2_x(eta,psi).^2+rv1_t2_y(eta,psi).^2+rv1_t2_z(eta,psi).^2).^(3/2);
f_xy_v2_1=@(eta,psi) 2*Av2_t1*(rv2_t1_x(eta,psi).*rv2_t1_y(eta,psi))./(rv2_t1_x(eta,psi).^2+rv2_t1_y(eta,psi).^2+rv2_t1_z(eta,psi).^2).^(3/2);
f_xy_v2_2=@(eta,psi) 2*Av2_t2*(rv2_t2_x(eta,psi).*rv2_t2_y(eta,psi))./(rv2_t2_x(eta,psi).^2+rv2_t2_y(eta,psi).^2+rv2_t2_z(eta,psi).^2).^(3/2);
f_xy_v3_1=@(eta,psi) 2*Av3_t1*(rv3_t1_x(eta,psi).*rv3_t1_y(eta,psi))./(rv3_t1_x(eta,psi).^2+rv3_t1_y(eta,psi).^2+rv3_t1_z(eta,psi).^2).^(3/2);
f_xy_v3_2=@(eta,psi) 2*Av3_t2*(rv3_t2_x(eta,psi).*rv3_t2_y(eta,psi))./(rv3_t2_x(eta,psi).^2+rv3_t2_y(eta,psi).^2+rv3_t2_z(eta,psi).^2).^(3/2);
f_xy_v4_1=@(eta,psi) 2*Av4_t1*(rv4_t1_x(eta,psi).*rv4_t1_y(eta,psi))./(rv4_t1_x(eta,psi).^2+rv4_t1_y(eta,psi).^2+rv4_t1_z(eta,psi).^2).^(3/2);
f_xy_v4_2=@(eta,psi) 2*Av4_t2*(rv4_t2_x(eta,psi).*rv4_t2_y(eta,psi))./(rv4_t2_x(eta,psi).^2+rv4_t2_y(eta,psi).^2+rv4_t2_z(eta,psi).^2).^(3/2);

% xz
f_xz_v1_1=@(eta,psi) 2*Av1_t1*(rv1_t1_x(eta,psi).*rv1_t1_z(eta,psi))./(rv1_t1_x(eta,psi).^2+rv1_t1_y(eta,psi).^2+rv1_t1_z(eta,psi).^2).^(3/2);
f_xz_v1_2=@(eta,psi) 2*Av1_t2*(rv1_t2_x(eta,psi).*rv1_t2_z(eta,psi))./(rv1_t2_x(eta,psi).^2+rv1_t2_y(eta,psi).^2+rv1_t2_z(eta,psi).^2).^(3/2);
f_xz_v2_1=@(eta,psi) 2*Av2_t1*(rv2_t1_x(eta,psi).*rv2_t1_z(eta,psi))./(rv2_t1_x(eta,psi).^2+rv2_t1_y(eta,psi).^2+rv2_t1_z(eta,psi).^2).^(3/2);
f_xz_v2_2=@(eta,psi) 2*Av2_t2*(rv2_t2_x(eta,psi).*rv2_t2_z(eta,psi))./(rv2_t2_x(eta,psi).^2+rv2_t2_y(eta,psi).^2+rv2_t2_z(eta,psi).^2).^(3/2);
f_xz_v3_1=@(eta,psi) 2*Av3_t1*(rv3_t1_x(eta,psi).*rv3_t1_z(eta,psi))./(rv3_t1_x(eta,psi).^2+rv3_t1_y(eta,psi).^2+rv3_t1_z(eta,psi).^2).^(3/2);
f_xz_v3_2=@(eta,psi) 2*Av3_t2*(rv3_t2_x(eta,psi).*rv3_t2_z(eta,psi))./(rv3_t2_x(eta,psi).^2+rv3_t2_y(eta,psi).^2+rv3_t2_z(eta,psi).^2).^(3/2);
f_xz_v4_1=@(eta,psi) 2*Av4_t1*(rv4_t1_x(eta,psi).*rv4_t1_z(eta,psi))./(rv4_t1_x(eta,psi).^2+rv4_t1_y(eta,psi).^2+rv4_t1_z(eta,psi).^2).^(3/2);
f_xz_v4_2=@(eta,psi) 2*Av4_t2*(rv4_t2_x(eta,psi).*rv4_t2_z(eta,psi))./(rv4_t2_x(eta,psi).^2+rv4_t2_y(eta,psi).^2+rv4_t2_z(eta,psi).^2).^(3/2);

% yy
f_yy_v1_1=@(eta,psi) 2*Av1_t1*(rv1_t1_y(eta,psi).*rv1_t1_y(eta,psi))./(rv1_t1_x(eta,psi).^2+rv1_t1_y(eta,psi).^2+rv1_t1_z(eta,psi).^2).^(3/2);
f_yy_v1_2=@(eta,psi) 2*Av1_t2*(rv1_t2_y(eta,psi).*rv1_t2_y(eta,psi))./(rv1_t2_x(eta,psi).^2+rv1_t2_y(eta,psi).^2+rv1_t2_z(eta,psi).^2).^(3/2);
f_yy_v2_1=@(eta,psi) 2*Av2_t1*(rv2_t1_y(eta,psi).*rv2_t1_y(eta,psi))./(rv2_t1_x(eta,psi).^2+rv2_t1_y(eta,psi).^2+rv2_t1_z(eta,psi).^2).^(3/2);
f_yy_v2_2=@(eta,psi) 2*Av2_t2*(rv2_t2_y(eta,psi).*rv2_t2_y(eta,psi))./(rv2_t2_x(eta,psi).^2+rv2_t2_y(eta,psi).^2+rv2_t2_z(eta,psi).^2).^(3/2);
f_yy_v3_1=@(eta,psi) 2*Av3_t1*(rv3_t1_y(eta,psi).*rv3_t1_y(eta,psi))./(rv3_t1_x(eta,psi).^2+rv3_t1_y(eta,psi).^2+rv3_t1_z(eta,psi).^2).^(3/2);
f_yy_v3_2=@(eta,psi) 2*Av3_t2*(rv3_t2_y(eta,psi).*rv3_t2_y(eta,psi))./(rv3_t2_x(eta,psi).^2+rv3_t2_y(eta,psi).^2+rv3_t2_z(eta,psi).^2).^(3/2);
f_yy_v4_1=@(eta,psi) 2*Av4_t1*(rv4_t1_y(eta,psi).*rv4_t1_y(eta,psi))./(rv4_t1_x(eta,psi).^2+rv4_t1_y(eta,psi).^2+rv4_t1_z(eta,psi).^2).^(3/2);
f_yy_v4_2=@(eta,psi) 2*Av4_t2*(rv4_t2_y(eta,psi).*rv4_t2_y(eta,psi))./(rv4_t2_x(eta,psi).^2+rv4_t2_y(eta,psi).^2+rv4_t2_z(eta,psi).^2).^(3/2);

% zy
f_zy_v1_1=@(eta,psi) 2*Av1_t1*(rv1_t1_z(eta,psi).*rv1_t1_y(eta,psi))./(rv1_t1_x(eta,psi).^2+rv1_t1_y(eta,psi).^2+rv1_t1_z(eta,psi).^2).^(3/2);
f_zy_v1_2=@(eta,psi) 2*Av1_t2*(rv1_t2_z(eta,psi).*rv1_t2_y(eta,psi))./(rv1_t2_x(eta,psi).^2+rv1_t2_y(eta,psi).^2+rv1_t2_z(eta,psi).^2).^(3/2);
f_zy_v2_1=@(eta,psi) 2*Av2_t1*(rv2_t1_z(eta,psi).*rv2_t1_y(eta,psi))./(rv2_t1_x(eta,psi).^2+rv2_t1_y(eta,psi).^2+rv2_t1_z(eta,psi).^2).^(3/2);
f_zy_v2_2=@(eta,psi) 2*Av2_t2*(rv2_t2_z(eta,psi).*rv2_t2_y(eta,psi))./(rv2_t2_x(eta,psi).^2+rv2_t2_y(eta,psi).^2+rv2_t2_z(eta,psi).^2).^(3/2);
f_zy_v3_1=@(eta,psi) 2*Av3_t1*(rv3_t1_z(eta,psi).*rv3_t1_y(eta,psi))./(rv3_t1_x(eta,psi).^2+rv3_t1_y(eta,psi).^2+rv3_t1_z(eta,psi).^2).^(3/2);
f_zy_v3_2=@(eta,psi) 2*Av3_t2*(rv3_t2_z(eta,psi).*rv3_t2_y(eta,psi))./(rv3_t2_x(eta,psi).^2+rv3_t2_y(eta,psi).^2+rv3_t2_z(eta,psi).^2).^(3/2);
f_zy_v4_1=@(eta,psi) 2*Av4_t1*(rv4_t1_z(eta,psi).*rv4_t1_y(eta,psi))./(rv4_t1_x(eta,psi).^2+rv4_t1_y(eta,psi).^2+rv4_t1_z(eta,psi).^2).^(3/2);
f_zy_v4_2=@(eta,psi) 2*Av4_t2*(rv4_t2_z(eta,psi).*rv4_t2_y(eta,psi))./(rv4_t2_x(eta,psi).^2+rv4_t2_y(eta,psi).^2+rv4_t2_z(eta,psi).^2).^(3/2);

% zz
f_zz_v1_1=@(eta,psi) 2*Av1_t1*(rv1_t1_z(eta,psi).*rv1_t1_z(eta,psi))./(rv1_t1_x(eta,psi).^2+rv1_t1_y(eta,psi).^2+rv1_t1_z(eta,psi).^2).^(3/2);
f_zz_v1_2=@(eta,psi) 2*Av1_t2*(rv1_t2_z(eta,psi).*rv1_t2_z(eta,psi))./(rv1_t2_x(eta,psi).^2+rv1_t2_y(eta,psi).^2+rv1_t2_z(eta,psi).^2).^(3/2);
f_zz_v2_1=@(eta,psi) 2*Av2_t1*(rv2_t1_z(eta,psi).*rv2_t1_z(eta,psi))./(rv2_t1_x(eta,psi).^2+rv2_t1_y(eta,psi).^2+rv2_t1_z(eta,psi).^2).^(3/2);
f_zz_v2_2=@(eta,psi) 2*Av2_t2*(rv2_t2_z(eta,psi).*rv2_t2_z(eta,psi))./(rv2_t2_x(eta,psi).^2+rv2_t2_y(eta,psi).^2+rv2_t2_z(eta,psi).^2).^(3/2);
f_zz_v3_1=@(eta,psi) 2*Av3_t1*(rv3_t1_z(eta,psi).*rv3_t1_z(eta,psi))./(rv3_t1_x(eta,psi).^2+rv3_t1_y(eta,psi).^2+rv3_t1_z(eta,psi).^2).^(3/2);
f_zz_v3_2=@(eta,psi) 2*Av3_t2*(rv3_t2_z(eta,psi).*rv3_t2_z(eta,psi))./(rv3_t2_x(eta,psi).^2+rv3_t2_y(eta,psi).^2+rv3_t2_z(eta,psi).^2).^(3/2);
f_zz_v4_1=@(eta,psi) 2*Av4_t1*(rv4_t1_z(eta,psi).*rv4_t1_z(eta,psi))./(rv4_t1_x(eta,psi).^2+rv4_t1_y(eta,psi).^2+rv4_t1_z(eta,psi).^2).^(3/2);
f_zz_v4_2=@(eta,psi) 2*Av4_t2*(rv4_t2_z(eta,psi).*rv4_t2_z(eta,psi))./(rv4_t2_x(eta,psi).^2+rv4_t2_y(eta,psi).^2+rv4_t2_z(eta,psi).^2).^(3/2);


% perform numerical integration
ymax = @(x) 1 - x;

Ires_1r_v1=quad2d(f_1r_v1_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
           quad2d(f_1r_v1_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_1r_v2=quad2d(f_1r_v2_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
           quad2d(f_1r_v2_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_1r_v3=quad2d(f_1r_v3_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
           quad2d(f_1r_v3_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_1r_v4=quad2d(f_1r_v4_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
           quad2d(f_1r_v4_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);

Ires_ij_v1(1)=quad2d(f_xx_v1_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
              quad2d(f_xx_v1_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_v2(1)=quad2d(f_xx_v2_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
              quad2d(f_xx_v2_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_v3(1)=quad2d(f_xx_v3_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
              quad2d(f_xx_v3_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_v4(1)=quad2d(f_xx_v4_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
              quad2d(f_xx_v4_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);

Ires_ij_v1(2)=quad2d(f_xy_v1_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
              quad2d(f_xy_v1_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_v2(2)=quad2d(f_xy_v2_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
              quad2d(f_xy_v2_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_v3(2)=quad2d(f_xy_v3_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
              quad2d(f_xy_v3_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_v4(2)=quad2d(f_xy_v4_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
              quad2d(f_xy_v4_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);

Ires_ij_v1(3)=quad2d(f_xz_v1_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
              quad2d(f_xz_v1_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_v2(3)=quad2d(f_xz_v2_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
              quad2d(f_xz_v2_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_v3(3)=quad2d(f_xz_v3_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
              quad2d(f_xz_v3_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_v4(3)=quad2d(f_xz_v4_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
              quad2d(f_xz_v4_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);

Ires_ij_v1(5)=quad2d(f_yy_v1_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
              quad2d(f_yy_v1_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_v2(5)=quad2d(f_yy_v2_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
              quad2d(f_yy_v2_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_v3(5)=quad2d(f_yy_v3_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
              quad2d(f_yy_v3_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_v4(5)=quad2d(f_yy_v4_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
              quad2d(f_yy_v4_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);

Ires_ij_v1(8)=quad2d(f_zy_v1_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
              quad2d(f_zy_v1_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_v2(8)=quad2d(f_zy_v2_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
              quad2d(f_zy_v2_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_v3(8)=quad2d(f_zy_v3_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
              quad2d(f_zy_v3_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_v4(8)=quad2d(f_zy_v4_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
              quad2d(f_zy_v4_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);

Ires_ij_v1(9)=quad2d(f_zz_v1_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
              quad2d(f_zz_v1_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_v2(9)=quad2d(f_zz_v2_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
              quad2d(f_zz_v2_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_v3(9)=quad2d(f_zz_v3_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
              quad2d(f_zz_v3_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);
Ires_ij_v4(9)=quad2d(f_zz_v4_1,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true)+...
              quad2d(f_zz_v4_2,0,1,0,ymax,'AbsTol',abs_tol,'Singular',true);

Ires_ij_v1(4)=Ires_ij_v1(2);
Ires_ij_v1(6)=Ires_ij_v1(8);
Ires_ij_v1(7)=Ires_ij_v1(3);

Ires_ij_v2(4)=Ires_ij_v2(2);
Ires_ij_v2(6)=Ires_ij_v2(8);
Ires_ij_v2(7)=Ires_ij_v2(3);

Ires_ij_v3(4)=Ires_ij_v3(2);
Ires_ij_v3(6)=Ires_ij_v3(8);
Ires_ij_v3(7)=Ires_ij_v3(3);

Ires_ij_v4(4)=Ires_ij_v4(2);
Ires_ij_v4(6)=Ires_ij_v4(8);
Ires_ij_v4(7)=Ires_ij_v4(3);
end