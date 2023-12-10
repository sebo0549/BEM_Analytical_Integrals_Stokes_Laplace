function M_quads = Calculate_Midpoints_Quads(Connectivity,Points)
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

% edge vectors of all quads
x1=Points(Connectivity(:,1),:);
x2=Points(Connectivity(:,2),:);
x3=Points(Connectivity(:,3),:);
x4=Points(Connectivity(:,4),:);

% calculate diagonals
d42=x4-x2;
d31=x3-x1;

% calc orientation regarding diagonal d42 
n11=Cross3D(d42,x3-x2);
n12=Cross3D(d42,x1-x2);

phi_normals=acosd(dot(n11,n12,2)./(sqrt(sum(n11.^2,2)).*sqrt(sum(n12.^2,2))));
ind_concave=phi_normals<90;

M_quads=x2+0.5*d42;
M_quads(ind_concave,:)=x1(ind_concave,:)+0.5*d31(ind_concave,:);
end

function N=Cross3D(a,b)
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

%vectorized version of cross product
N = [a(:,2).*b(:,3) - a(:,3).*b(:,2) a(:,3).*b(:,1) - a(:,1).*b(:,3) a(:,1).*b(:,2) - a(:,2).*b(:,1)]; 
end