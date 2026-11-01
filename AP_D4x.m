function D4x    =   AP_D4x(N,h,BC)
% Author: Abhiroop Paul
% Date: 30th July 2026
%
% Title: 4th Order Spatial Difference Operator Matrix
% 
% Description:
% Function to calculate the 4th order spatial difference operator as a
% matrix. The output matrix encodes the boundary condition.
% Input argument needs the size of vector the matrix is to be acted on "N",
% the spatial sampling step "h" and the boundary condition ID "BC" where:
% BC = 0 => clamped
% BC = 1 => simply suppported
% BC = 2 => free
% =======================================================================

% Spatial 4th Order Difference with BCs encoded
D4x     =   ((1/(h^4))*toeplitz([6;-4;1;zeros(N+1-3,1)]));
D4x(1,:)=   zeros(N+1,1); 
D4x(N+1,:)= zeros(N+1,1);

if floor(mod(BC,3)) == 0
    % Clamped
    D4x(2,:)=   zeros(N+1,1);
    D4x(N,:)=   zeros(N+1,1);
elseif floor(mod(BC,3)) == 1
    % Simply Supported
    D4x(2,:)=   (1/(h^4))*[0;5;-4;1;zeros(N+1-4,1)];
    D4x(N,:)=   (1/(h^4))*[zeros(N+1-4,1);1;-4;5;0];
else
    % Free
    D4x(1,:)    = (1/(h^4))*[2; -4; 2; zeros(N+1-3,1)];
    D4x(2,:)    = (1/(h^4))*[-2; 5; -4; 1; zeros(N+1-4,1)];
    D4x(N,:)    = (1/(h^4))*[zeros(N+1-4,1); 1; -4; 5; -2];
    D4x(N+1,:)  = (1/(h^4))*[zeros(N+1-3,1); 2; -4; 2];
end

end
