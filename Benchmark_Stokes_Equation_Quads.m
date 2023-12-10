function Benchmark_Stokes_Equation_Quads
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

close all
addpath(genpath('Functions'))
addpath(genpath('Inputfiles'))

% number of runs for the benchmark
Nbenchmark=1000;

% import quad mesh
fprintf('Load quad mesh\n')
[P,Con] = Import_Comsol_Mesh_File_Quads('test_geometry_quads_100t_elements.mphtxt');

% calculate midpoints of all quads
M_Quads = Calculate_Midpoints_Quads(Con,P);

% perform benchmark
fprintf('Calculate integrals:\n')
tn_midpoints=zeros(Nbenchmark,1);
tn_vertices=zeros(Nbenchmark,1);
for n=1:Nbenchmark
    tic
    Calculate_Integrals_Stokes_Flow_Quads_Midpoints(P,M_Quads,Con);
    tn_midpoints(n)=toc;

    tic
    Calculate_Integrals_Stokes_Flow_Quads_Vertices(P,Con);
    tn_vertices(n)=toc;

    % print time needed
    fprintf('   %1.0f. Benchmark run\n',n)
    fprintf('   Time needed using inner points: t_%1.0d=%1.3f ms\n',n,tn_midpoints(n)*1e3)
    fprintf('   Time needed using vertices:     t_%1.0d=%1.3f ms\n\n',n,tn_vertices(n)*1e3)
end

% calculate mean calculation time and the standard deviation
t_mean_Midpoints=mean(tn_midpoints);
t_std_Midpoints=std(tn_midpoints);

t_mean_Vertices=mean(tn_vertices);
t_std_Vertices=std(tn_vertices);

% calc number of elements per second
num_quads_per_sec_midpoints=length(Con(:,1))/t_mean_Midpoints;
num_quads_per_sec_vertices=length(Con(:,1))/t_mean_Vertices;

% print results
fprintf('Results:\n')
fprintf('Number of elements per run: %1.0d\n',length(Con(:,1)))
fprintf('Total number of elements: %1.0d\n',Nbenchmark*length(Con(:,1)))
fprintf('Total time needed using inner points: t=%1.3f s\n',sum(tn_midpoints))
fprintf('Total time needed using vertices:     t=%1.3f s\n',sum(tn_vertices))
fprintf('Average time needed per run using inner points: t=(%1.3f+-%1.3f) ms\n',t_mean_Midpoints*1e3,t_std_Midpoints*1e3)
fprintf('Average time needed per run using vertices:     t=(%1.3f+-%1.3f) ms\n',t_mean_Vertices*1e3,t_std_Vertices*1e3)
fprintf('Average number of elements per second using inner points: %1.3d\n',num_quads_per_sec_midpoints)
fprintf('Average number of elements per second using vertices:     %1.3d\n',num_quads_per_sec_vertices)

% plot quad mesh
figure()
hold on
patch('faces',Con,'vertices',P,'facecolor','none','FaceAlpha',0.1)
view(3)
grid on
axis equal
end







