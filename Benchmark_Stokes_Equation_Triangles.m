function Benchmark_Stokes_Equation_Triangles
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

close all
addpath(genpath('Functions'))
addpath(genpath('Inputfiles'))

% number of runs for the benchmark
Nbenchmark=1000;

% import triangulation
fprintf('Load triangulation\n')
TR=stlread('test_geometry_triangles_500t_elements.stl');

% extract points and connectivity list
P=TR.Points;
Con=TR.ConnectivityList;

% calculate incenters and centroids of all triangles
[M_Incenter,M_Centroid] = Calculate_Midpoints_Triangles(Con,P);

% perform benchmark
fprintf('Calculate integrals:\n')
tn_Incenter=zeros(Nbenchmark,1);
tn_Centroid=zeros(Nbenchmark,1);
tn_Vertices=zeros(Nbenchmark,1);
for n=1:Nbenchmark
    tic
    Calculate_Integrals_Stokes_Flow_Triangles_Incenter(P,M_Incenter,Con);
    tn_Incenter(n)=toc;

    tic
    Calculate_Integrals_Stokes_Flow_Triangles_Centroid(P,M_Centroid,Con);
    tn_Centroid(n)=toc;

    tic
    Calculate_Integrals_Stokes_Flow_Triangles_Vertices(P,Con);
    tn_Vertices(n)=toc;
    
    % print time needed
    fprintf('   %1.0f. Benchmark run\n',n)
    fprintf('   Time needed using incenters t_%1.0d=%1.3f ms\n',n,tn_Incenter(n)*1e3)
    fprintf('   Time needed using centroids t_%1.0d=%1.3f ms\n',n,tn_Centroid(n)*1e3)
    fprintf('   Time needed using vertices  t_%1.0d=%1.3f ms\n\n',n,tn_Vertices(n)*1e3)
end

% calculate mean calculation time and the standard deviation
t_mean_Incenter=mean(tn_Incenter);
t_std_Incenter=std(tn_Incenter);

t_mean_Centroid=mean(tn_Centroid);
t_std_Centroid=std(tn_Centroid);

t_mean_Vertices=mean(tn_Vertices);
t_std_Vertices=std(tn_Vertices);

% print number of elements per second
num_tris_per_sec_Incenter=length(Con(:,1))/t_mean_Incenter;
num_tris_per_sec_Centroid=length(Con(:,1))/t_mean_Centroid;
num_tris_per_sec_Vertices=length(Con(:,1))/t_mean_Vertices;

% print results
fprintf('Results:\n')
fprintf('Number of elements per run: %1.0d\n',length(Con(:,1)))
fprintf('Total number of elements:   %1.0d\n',Nbenchmark*length(Con(:,1)))
fprintf('Total time needed using incenters: t=%1.3f s\n',sum(tn_Incenter))
fprintf('Total time needed using centroids: t=%1.3f s\n',sum(tn_Centroid))
fprintf('Total time needed using vertices:  t=%1.3f s\n',sum(tn_Vertices))
fprintf('Average time needed per run using incenters: t=(%1.3f+-%1.3f) ms\n',t_mean_Incenter*1e3,t_std_Incenter*1e3)
fprintf('Average time needed per run using centroids: t=(%1.3f+-%1.3f) ms\n',t_mean_Centroid*1e3,t_std_Centroid*1e3)
fprintf('Average time needed per run using vertices:  t=(%1.3f+-%1.3f) ms\n',t_mean_Vertices*1e3,t_std_Vertices*1e3)
fprintf('Average number of elements per second using incenters: %1.3d\n',num_tris_per_sec_Incenter)
fprintf('Average number of elements per second using centroids: %1.3d\n',num_tris_per_sec_Centroid)
fprintf('Average number of elements per second using vertices:  %1.3d\n',num_tris_per_sec_Vertices)

% plot triangulation
figure()
set(gcf,'units','normalized','outerposition',[0 0 1 1])
trisurf(Con,P(:,1),P(:,2),P(:,3),'FaceColor',[0.5 0.5 0.5],'EdgeAlpha',0.25)
grid on
axis equal
drawnow
end







