function Comparison_Numerical_Integrals_Triangles
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

close all
format long

addpath(genpath('Functions'))
addpath(genpath('Inputfiles'))

abs_tol=1e-12; % tolerance for numerical solution

%% Define triangulation using a connectivity and a point list
Con=[1,2,3];
P  =[-0.12,0,0; 0,0.1,0.23077; 0,0.077,0];

%% Calculate incenters and centroids of all triangles
[M_Incenter,M_Centroid] = Calculate_Midpoints_Triangles(Con,P);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analytical calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the singular integrals for the Stokes flow using analytical equations
[Ires_ij_v1,Ires_1r_v1,...
    Ires_ij_v2,Ires_1r_v2,...
    Ires_ij_v3,Ires_1r_v3]=Calculate_Integrals_Stokes_Flow_Triangles_Vertices(P,Con);
[Ires_ij_Incenter,Ires_1r_Incenter]=Calculate_Integrals_Stokes_Flow_Triangles_Incenter(P,M_Incenter,Con);
[Ires_ij_Centroid,Ires_1r_Centroid]=Calculate_Integrals_Stokes_Flow_Triangles_Centroid(P,M_Centroid,Con);

% Calculate the singular integrals for the Laplace equation using analytical equations
[I_1r_v1_analytical_laplace,I_1r_v2_analytical_laplace,I_1r_v3_analytical_laplace]=Calculate_Integrals_Laplace_Triangles_Vertices(P,Con);
I_1r_Incenter_analytical_laplace=Calculate_Integrals_Laplace_Triangles_Incenter(P,M_Incenter,Con);
I_1r_Centroid_analytical_laplace=Calculate_Integrals_Laplace_Triangles_Centroid(P,M_Centroid,Con);

% Reshape solutions for stokes flow
% incenter - reshape analytical integrals into matrix form 
[I_analytical_Incenter,F_analytical_Incenter,H_analytical_Incenter]=...
    Reshape_Result_Analytical(Ires_ij_Incenter,Ires_1r_Incenter);

% centroids - reshape analytical integrals into matrix form 
[I_analytical_Centroid,F_analytical_Centroid,H_analytical_Centroid]=...
    Reshape_Result_Analytical(Ires_ij_Centroid,Ires_1r_Centroid);

% 1. vertex - reshape analytical integrals into matrix form 
[I_analytical_v1,F_analytical_v1,H_analytical_v1]=...
    Reshape_Result_Analytical(Ires_ij_v1,Ires_1r_v1);

% 2. vertex - reshape analytical integrals into matrix form 
[I_analytical_v2,F_analytical_v2,H_analytical_v2]=...
    Reshape_Result_Analytical(Ires_ij_v2,Ires_1r_v2);

% 3. vertex - reshape analytical integrals into matrix form 
[I_analytical_v3,F_analytical_v3,H_analytical_v3]=...
    Reshape_Result_Analytical(Ires_ij_v3,Ires_1r_v3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Numerical calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the singular integrals for the Stokes flow and the Laplace equation using a numerical approximation
[Ires_ij_v1_num,Ires_1r_v1_num,...
    Ires_ij_v2_num,Ires_1r_v2_num,...
    Ires_ij_v3_num,Ires_1r_v3_num]=Calculate_Integrals_Numerical_Triangles_Vertices(P,abs_tol);
[Ires_ij_Incenter_num,Ires_1r_Incenter_num]=Calculate_Integrals_Numerical_Triangles_Midpoints(P,M_Incenter,abs_tol);
[Ires_ij_Centroid_num,Ires_1r_Centroid_num]=Calculate_Integrals_Numerical_Triangles_Midpoints(P,M_Centroid,abs_tol);

% Reshape solutions for stokes flow
% incenter - reshape numerical integrals into matrix form 
[I_num_Incenter,F_num_Incenter,H_num_Incenter]=...
    Reshape_Result_Numerical(Ires_ij_Incenter_num,Ires_1r_Incenter_num);

% centroids - reshape numerical integrals into matrix form 
[I_num_Centroid,F_num_Centroid,H_num_Centroid]=...
    Reshape_Result_Numerical(Ires_ij_Centroid_num,Ires_1r_Centroid_num);

% 1. vertex - reshape numerical integrals into matrix form 
[I_num_v1,F_num_v1,H_num_v1]=...
    Reshape_Result_Numerical(Ires_ij_v1_num,Ires_1r_v1_num);

% 2. vertex - reshape numerical integrals into matrix form 
[I_num_v2,F_num_v2,H_num_v2]=...
    Reshape_Result_Numerical(Ires_ij_v2_num,Ires_1r_v2_num);

% 3. vertex - reshape numerical integrals into matrix form 
[I_num_v3,F_num_v3,H_num_v3]=...
    Reshape_Result_Numerical(Ires_ij_v3_num,Ires_1r_v3_num);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% print results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Stokes equation:')
disp('H - Incenter')
H_num_Incenter
H_analytical_Incenter

disp('H - Centroid')
H_num_Centroid
H_analytical_Centroid

disp('H - v1')
H_num_v1
H_analytical_v1

disp('H - v2')
H_num_v2
H_analytical_v2

disp('H - v3')
H_num_v3
H_analytical_v3

disp('F - Incenter')
F_num_Incenter
F_analytical_Incenter

disp('F - Centroid')
F_num_Centroid
F_analytical_Centroid

disp('F - v1')
F_num_v1
F_analytical_v1

disp('F - v2')
F_num_v2
F_analytical_v2

disp('F - v3')
F_num_v3
F_analytical_v3

disp('I - Incenter')
I_num_Incenter
I_analytical_Incenter   

disp('I - Centroid')
I_num_Centroid
I_analytical_Centroid   

disp('I - v1')
I_num_v1
I_analytical_v1 

disp('I - v2')
I_num_v2
I_analytical_v2 

disp('I - v3')
I_num_v3
I_analytical_v3 

disp('Laplace equation:')
Ires_1r_Incenter_num
I_1r_Incenter_analytical_laplace

Ires_1r_Centroid_num
I_1r_Centroid_analytical_laplace

Ires_1r_v1_num
I_1r_v1_analytical_laplace

Ires_1r_v2_num
I_1r_v2_analytical_laplace

Ires_1r_v3_num
I_1r_v3_analytical_laplace


%% plot triangle together with the midpoints
figure()
hold on
plot3(P(:,1),P(:,2),P(:,3),'.k','MarkerSize',15)
plot3(M_Incenter(1),M_Incenter(2),M_Incenter(3),'.g','MarkerSize',12)
plot3(M_Centroid(1),M_Centroid(2),M_Centroid(3),'.r','MarkerSize',12)
trisurf(Con,P(:,1),P(:,2),P(:,3),'facecolor','blue','FaceAlpha',0.1)
grid on
view(0,0)
axis equal
end









