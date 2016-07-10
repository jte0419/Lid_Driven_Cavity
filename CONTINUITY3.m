function [RHSc] = CONTINUITY3(N,uu,vv,zeta_x,eta_y,J,dZeta,dEta,D1c,D2c)

% INPUTS
% - N      : Number of grid nodes
% - uu     : U velocity field [NxN]
% - vv     : V velocity field [NxN]
% - zeta_x : Metric of the transformation in zeta direction [NxN]
% - eta_y  : Metric of the transformation in eta direction [NxN]
% - J      : Jacobian of the transformation [NxN]
% - dZeta  : Step size in zeta direction of computational space
% - dEta   : Step size in eta direction of computational space
% - D1c    : Dissipation in zeta direction [NxN-1]
% - D2c    : Dissipation in eta direction [N-1xN]
% 
% OUTPUTS
% - RHSc : Right-hand-side matrix for continuity [N-2xN-2]

% Initialize matrices
dD1c_dZ = zeros(N-2,N-2);
dD2c_dE = zeros(N-2,N-2);
dE1c_dZ = zeros(N-2,N-2);
dE2c_dE = zeros(N-2,N-2);
Rc_n    = zeros(N-2,N-2);

% Inviscid Flux Terms
% - [25x25] matrix
E1c = (uu.*zeta_x)./J;
E2c = (vv.*eta_y)./J;

% Dissipation central finite difference
% - D1c is [25x24]
% - D2c is [24x25]
for row = 1:1:N-2
    for col = 1:1:N-2
        dD1c_dZ(row,col) = (D1c(row+1,col+1)-D1c(row+1,col))/dZeta;
    end
end
for row = 1:1:N-2
    for col = 1:1:N-2
        dD2c_dE(row,col) = (D2c(row+1,col+1)-D2c(row,col+1))/dEta;
    end
end
assignin('base','dD1c_dZ',dD1c_dZ);
assignin('base','dD2c_dE',dD2c_dE);

% Calculate finite differences and obtain RHS of continuity
% - Don't need to evaluate on the boundaries so: [N-2xN-2]
% - Jacobian is full [NxN] so start at i+1 and j+1
% - E1c and E2c are [25x25] matrices
for row = 1:1:N-2
    for col = 1:1:N-2
        dE1c_dZ(row,col) = (E1c(row+1,col+2)-E1c(row+1,col))/(2*dZeta);
        dE2c_dE(row,col) = (E2c(row+2,col+1)-E2c(row,col+1))/(2*dEta);
    end
end

for row = 1:1:N-2
    for col = 1:1:N-2
        Rc_n(row,col) = J(row+1,col+1)*(-dE1c_dZ(row,col)-dE2c_dE(row,col)-...
                                dD1c_dZ(row,col)-dD2c_dE(row,col));
    end
end

RHSc = Rc_n;

assignin('base','E1c',E1c);
assignin('base','E2c',E2c);
assignin('base','dE1c_dZ',dE1c_dZ);
assignin('base','dE2c_dE',dE2c_dE);
