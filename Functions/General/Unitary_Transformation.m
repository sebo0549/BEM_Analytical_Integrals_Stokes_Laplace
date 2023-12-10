function [e1s,e2s]=Unitary_Transformation(b1s,b2s,L)
% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 07-12-2023
% License: MIT License

% calculate the new basis vectors
e1s=b1s./L;
e2s=b2s./sqrt(sum(b2s.^2,2));
end