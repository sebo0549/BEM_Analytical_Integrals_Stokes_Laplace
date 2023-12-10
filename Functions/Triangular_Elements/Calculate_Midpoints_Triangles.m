function [M_IC,M_Cent] = Calculate_Midpoints_Triangles(Connectivity,Points)
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

% calculate centroid coordinates of all triangles
M_Cent=(1/3)*(Points(Connectivity(:,1),:)+...
              Points(Connectivity(:,2),:)+...
              Points(Connectivity(:,3),:));

%calculate incenter coordinates of all triangles
d12=Points(Connectivity(:,2),:)-Points(Connectivity(:,1),:);
d13=Points(Connectivity(:,3),:)-Points(Connectivity(:,1),:);
d32=Points(Connectivity(:,3),:)-Points(Connectivity(:,2),:);

L12=sqrt(sum(d12.^2,2));
L13=sqrt(sum(d13.^2,2));
L32=sqrt(sum(d32.^2,2));

M_IC=(Points(Connectivity(:,3),:).*L12+...
      Points(Connectivity(:,2),:).*L13+...
      Points(Connectivity(:,1),:).*L32)./(L12+L13+L32);
end