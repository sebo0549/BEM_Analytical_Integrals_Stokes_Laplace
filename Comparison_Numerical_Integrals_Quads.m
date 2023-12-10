function Comparison_Numerical_Integrals_Quads
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

close all
format long

addpath(genpath('Functions'))
addpath(genpath('Inputfiles'))

abs_tol=1e-12; % tolerance for numerical solution

%% Define quad mesh using a connectivity and a point list
P=[-0.5,-1.7,0; 0.25,-2.5,0; 1,-0.75,0; 0,1,0];
Con=[1,2,3,4];

%% Calculate midpoints of all quads
M_Quads = Calculate_Midpoints_Quads(Con,P);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analytical calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the singular integrals for the Stokes flow using analytical equations
[Ires_ij_v1,Ires_1r_v1,...
    Ires_ij_v2,Ires_1r_v2,...
    Ires_ij_v3,Ires_1r_v3,...
    Ires_ij_v4,Ires_1r_v4]=Calculate_Integrals_Stokes_Flow_Quads_Vertices(P,Con);
[Ires_ij_Midpoints,Ires_1r_Midpoints]=Calculate_Integrals_Stokes_Flow_Quads_Midpoints(P,M_Quads,Con);

% Calculate the singular integrals for the Laplace equation using analytical equations
[I_1r_v1_analytical_laplace,I_1r_v2_analytical_laplace,...
    I_1r_v3_analytical_laplace,I_1r_v4_analytical_laplace]=Calculate_Integrals_Laplace_Quads_Vertices(P,Con);
I_1r_Midpoints_analytical_laplace=Calculate_Integrals_Laplace_Quads_Midpoints(P,M_Quads,Con);

% Reshape solutions for stokes flow
% midpoints - reshape analytical integrals into matrix form 
[I_analytical_Midpoints,F_analytical_Midpoints,H_analytical_Midpoints]=...
    Reshape_Result_Analytical(Ires_ij_Midpoints,Ires_1r_Midpoints);

% 1. vertex - reshape analytical integrals into matrix form 
[I_analytical_v1,F_analytical_v1,H_analytical_v1]=...
    Reshape_Result_Analytical(Ires_ij_v1,Ires_1r_v1);

% 2. vertex - reshape analytical integrals into matrix form 
[I_analytical_v2,F_analytical_v2,H_analytical_v2]=...
    Reshape_Result_Analytical(Ires_ij_v2,Ires_1r_v2);

% 3. vertex - reshape analytical integrals into matrix form 
[I_analytical_v3,F_analytical_v3,H_analytical_v3]=...
    Reshape_Result_Analytical(Ires_ij_v3,Ires_1r_v3);

% 4. vertex - reshape analytical integrals into matrix form 
[I_analytical_v4,F_analytical_v4,H_analytical_v4]=...
    Reshape_Result_Analytical(Ires_ij_v4,Ires_1r_v4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Numerical calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the singular integrals for the Stokes flow using a numerical approximation
[Ires_ij_Midpoints_num,Ires_1r_Midpoints_num]=Calculate_Integrals_Numerical_Quads_Midpoints(P,M_Quads,abs_tol);
[Ires_ij_v1_num,Ires_1r_v1_num,...
    Ires_ij_v2_num,Ires_1r_v2_num,...
    Ires_ij_v3_num,Ires_1r_v3_num,...
    Ires_ij_v4_num,Ires_1r_v4_num]=Calculate_Integrals_Numerical_Quads_Vertices(P,abs_tol);

% Reshape solutions for stokes flow
% midpoints - reshape numerical integrals into matrix form 
[I_num_Midpoints,F_num_Midpoints,H_num_Midpoints]=...
    Reshape_Result_Numerical(Ires_ij_Midpoints_num,Ires_1r_Midpoints_num);

% 1. vertex - reshape numerical integrals into matrix form 
[I_num_v1,F_num_v1,H_num_v1]=...
    Reshape_Result_Numerical(Ires_ij_v1_num,Ires_1r_v1_num);

% 2. vertex - reshape numerical integrals into matrix form 
[I_num_v2,F_num_v2,H_num_v2]=...
    Reshape_Result_Numerical(Ires_ij_v2_num,Ires_1r_v2_num);

% 3. vertex - reshape numerical integrals into matrix form 
[I_num_v3,F_num_v3,H_num_v3]=...
    Reshape_Result_Numerical(Ires_ij_v3_num,Ires_1r_v3_num);

% 4. vertex - reshape numerical integrals into matrix form 
[I_num_v4,F_num_v4,H_num_v4]=...
    Reshape_Result_Numerical(Ires_ij_v4_num,Ires_1r_v4_num);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% print results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Stokes equation:')
disp('H - Midpoints')
H_num_Midpoints
H_analytical_Midpoints

disp('H - v1')
H_num_v1
H_analytical_v1

disp('H - v2')
H_num_v2
H_analytical_v2

disp('H - v3')
H_num_v3
H_analytical_v3

disp('H - v4')
H_num_v4
H_analytical_v4

disp('F - Midpoints')
F_num_Midpoints
F_analytical_Midpoints

disp('F - v1')
F_num_v1
F_analytical_v1

disp('F - v2')
F_num_v2
F_analytical_v2

disp('F - v3')
F_num_v3
F_analytical_v3

disp('F - v4')
F_num_v4
F_analytical_v4


disp('I - Midpoints')
I_num_Midpoints
I_analytical_Midpoints  

disp('I - v1')
I_num_v1
I_analytical_v1

disp('I - v2')
I_num_v2
I_analytical_v2

disp('I - v3')
I_num_v3
I_analytical_v3

disp('I - v4')
I_num_v4
I_analytical_v4

disp('Laplace equation:')
Ires_1r_Midpoints_num
I_1r_Midpoints_analytical_laplace

Ires_1r_v1_num
I_1r_v1_analytical_laplace

Ires_1r_v2_num
I_1r_v2_analytical_laplace

Ires_1r_v3_num
I_1r_v3_analytical_laplace

Ires_1r_v4_num
I_1r_v4_analytical_laplace


%% plot quad together with the midpoint
figure()
hold on
plot3(P(:,1),P(:,2),P(:,3),'.k','MarkerSize',15)
plot3(M_Quads(:,1),M_Quads(:,2),M_Quads(:,3),'.r','MarkerSize',15)
patch('faces',Con,'vertices',P,'facecolor','blue','FaceAlpha',0.1)
view(2)
grid on
axis equal
end







