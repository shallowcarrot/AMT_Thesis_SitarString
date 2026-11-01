function D2x    =   AP_D2x(N,h,BC)
% Author: Abhiroop Paul
% Date: 30th July 2026
%
% Title: 2nd Order Spatial Difference Matrix
% 
% Description:
% Function to calculate the 2nd order spatial difference operator as a
% matrix. The output matrix encodes the boundary condition.
% Input argument needs the size of vector the matrix is to be acted on "N",
% the spatial sampling step "h" and the boundaray condition ID "BC" where:
% BC = 0 => clamped
% BC = 1 => simply suppported
% BC = 2 => free
% =======================================================================

% Spatial 2nd Order Difference with first and last rows = 0 to satisfy
% boundary conditions
D2x         =   ((1/(h^2))*toeplitz([-2;1;zeros(N+1-2,1)]));
D2x(1,:)    =   zeros(N+1,1);
D2x(N+1,:)  =   zeros(N+1,1);

end