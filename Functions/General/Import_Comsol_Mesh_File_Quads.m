function [Points,Connectivity] = Import_Comsol_Mesh_File_Quads(file_name)
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

mesh_file=textread(file_name,'%s','delimiter','\n');

% empty lines
rows_empty=find(strcmp(mesh_file,''));

% all lines including "# Elements"
rows_elements=find(strcmp(mesh_file,'# Elements'));

% index of first row of the point coordinates -> comsol v55
rows_xf_start=find(strcmp(mesh_file,'# Mesh point coordinates'))+1;

if(isempty(rows_xf_start))  % -> comsol v56
    rows_xf_start=find(strcmp(mesh_file,'# Mesh vertex coordinates'))+1;
end

% index of last row of the point coordinates
rows_xf_end=rows_empty(find(rows_empty>rows_xf_start,1,'first'));

% index of first row of the triangles
rows_quads_start=find(strcmp(mesh_file,'4 quad # type name'));

% index of the first row of the connectivity list
rows_quads_start_elem=rows_elements(find(rows_elements>rows_quads_start,1))+1;
% index of the last row of the connectivity list
rows_quads_end_elem=rows_empty(find(rows_empty>rows_quads_start_elem,1));

% read coordinates of the mesh vertices
for n=1:(rows_xf_end-rows_xf_start)
    Points(n,:)=sscanf(mesh_file{n+rows_xf_start-1}, '%f%f%f', 3)';
end

% read connectivity of the mesh vertices
for n=1:(rows_quads_end_elem-rows_quads_start_elem)
    Connectivity(n,:)=sscanf(mesh_file{n+rows_quads_start_elem-1}, '%d%d%d%d', 4)'+1;
end
Connectivity=[Connectivity(:,2) Connectivity(:,1) Connectivity(:,3) Connectivity(:,4)];
end

