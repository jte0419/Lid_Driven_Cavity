function [RHSu] = XMOMENTUM3(N,uu,vv,PP,zeta_x,eta_y,J,dZeta,dEta,Re,D1u,D2u)

% INPUTS
% - N      : Number of grid nodes
% - uu     : U velocity field [N,N]
% - vv     : V velocity field [N,N]
% - zeta_x : Metric of the transformation in zeta direction [N,N]
% - eta_y  : Metric of the transformation in eta direction [N,N]
% - J      : Jacobian of the transformation [N,N]
% - dZeta  : Step size in zeta direction of computational space
% - dEta   : Step size in eta direction of computational space
% - D1u    : Dissipation in zeta direction [N,N-1]
% - D2u    : Dissipation in eta direction [N-1,N]
% 
% OUTPUTS
% - RHSu : Right-hand-side matrix for x-momentum [N-2,N-2]

% Initialize matrices
dE1u_dZ   = zeros(N-2,N-2);
dE2u_dE   = zeros(N-1,N-2);
E1vu_hfnd = zeros(N,N-1);
E2vu_hfnd = zeros(N-1,N);
dE1vu_dZ  = zeros(N-2,N-2);
dE2vu_dE  = zeros(N-2,N-2);
dD1u_dZ   = zeros(N-2,N-2);
dD2u_dE   = zeros(N-2,N-2);
Ru_n      = zeros(N-2,N-2);

% Inviscid Terms
% - [N,N] matrices
E1u = ((uu.*uu.*zeta_x)+(PP.*zeta_x))./J;
E2u = (uu.*vv.*eta_y)./J;

% Inviscid Derivatives
% - E1u and E2u are [N,N] matrices
for row = 1:1:N-2
    for col = 1:1:N-2
        dE1u_dZ(row,col) = (E1u(row+1,col+2)-E1u(row+1,col))/(2*dZeta);
        dE2u_dE(row,col) = (E2u(row+2,col+1)-E2u(row,col+1))/(2*dEta);
    end
end

% Viscous Zeta Term (halfnode)
% [N,N-1]
for row = 1:1:N
    for col = 1:1:N-1
        zeta_x_hfnd  = (1/2)*(zeta_x(row,col+1)+zeta_x(row,col));
        g11_hfnd     = zeta_x_hfnd^2;
        du_dZ_hfnd   = (uu(row,col+1)-uu(row,col))/(dZeta);
        J_x_hfnd     = (1/2)*(J(row,col+1)+J(row,col));
        
        % Halfnode viscous zeta term
        E1vu_hfnd(row,col) = (g11_hfnd/(Re*J_x_hfnd))*...
                                (du_dZ_hfnd);
    end
end
assignin('base','E1vu_hfnd',E1vu_hfnd);

% Viscous Eta Term (halfnode)
% [N-1,N]
for row = 1:1:N-1
    for col = 1:1:N
        eta_y_hfnd   = (1/2)*(eta_y(row+1,col)+eta_y(row,col));
        g22_hfnd     = eta_y_hfnd^2;
        du_dE_hfnd   = (uu(row+1,col)-uu(row,col))/(dEta);
        J_y_hfnd     = (1/2)*(J(row+1,col)+J(row,col)); 

        % Halfnode viscous eta term
        E2vu_hfnd(row,col) = (g22_hfnd/(Re*J_y_hfnd))*...
                                (du_dE_hfnd);        
    end
end
assignin('base','E2vu_hfnd',E2vu_hfnd);

% Viscous Derivatives
% - E1vu_hfnd is [N,N-1] matrix
% - E2vu_hfnd is [N-1,N] matrix
for row = 1:1:N-2
    for col = 1:1:N-2
        dE1vu_dZ(row,col) = (E1vu_hfnd(row+1,col+1)-E1vu_hfnd(row+1,col))/dZeta;
        dE2vu_dE(row,col) = (E2vu_hfnd(row+1,col+1)-E2vu_hfnd(row,col+1))/dEta;
    end
end

% Dissipation central finite difference
% - D1u is [25x24]
% - D2u is [24x25]
for row = 1:1:N-2
    for col = 1:1:N-2
        dD1u_dZ(row,col) = (D1u(row+1,col+1)-D1u(row+1,col))/dZeta;
    end
end
for row = 1:1:N-2
    for col = 1:1:N-2
        dD2u_dE(row,col) = (D2u(row+1,col+1)-D2u(row,col+1))/dEta;
    end
end
assignin('base','dD1u_dZ',dD1u_dZ);
assignin('base','dD2u_dE',dD2u_dE);

for row = 1:1:N-2
    for col = 1:1:N-2
        Ru_n(row,col) = J(row+1,col+1)*(-dE1u_dZ(row,col)-dE2u_dE(row,col)+...
                            dE1vu_dZ(row,col)+dE2vu_dE(row,col)-...
                            dD1u_dZ(row,col)-dD2u_dE(row,col));
    end
end

RHSu = Ru_n;

assignin('base','E1u',E1u);
assignin('base','E2u',E2u);
assignin('base','dE1u_dZ',dE1u_dZ);
assignin('base','dE2u_dE',dE2u_dE);

