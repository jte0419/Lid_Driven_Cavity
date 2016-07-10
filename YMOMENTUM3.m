function [RHSv] = YMOMENTUM3(N,uu,vv,PP,zeta_x,eta_y,J,dZeta,dEta,Re,D1v,D2v)

% INPUTS
% - N      : Number of grid nodes
% - uu     : U velocity field [NxN]
% - vv     : V velocity field [NxN]
% - zeta_x : Metric of the transformation in zeta direction [NxN]
% - eta_y  : Metric of the transformation in eta direction [NxN]
% - J      : Jacobian of the transformation [NxN]
% - dZeta  : Step size in zeta direction of computational space
% - dEta   : Step size in eta direction of computational space
% - D1v    : Dissipation in zeta direction [NxN-1]
% - D2v    : Dissipation in eta direction [N-1xN]
% 
% OUTPUTS
% - RHSv : Right-hand-side matrix for y-momentum [N-2xN-2]

% Initialize matrices
dE1v_dZ   = zeros(N-2,N-2);
dE2v_dE   = zeros(N-2,N-2);
E1vv_hfnd = zeros(N,N-1);
E2vv_hfnd = zeros(N-1,N);
dD1v_dZ   = zeros(N-2,N-2);
dD2v_dE   = zeros(N-2,N-2);
dE1vv_dZ  = zeros(N-2,N-2);
dE2vv_dE  = zeros(N-2,N-2);
Rv_n      = zeros(N-2,N-2);

% Inviscid Terms
% - [25x25] matrices
E1v = (vv.*uu.*zeta_x)./J;
E2v = ((vv.*vv.*eta_y)+(PP.*eta_y))./J;

% Inviscid Derivatives
% - E1v and E2v are [25x25] matrices
for row = 1:1:N-2
    for col = 1:1:N-2
        dE1v_dZ(row,col) = (E1v(row+1,col+2)-E1v(row+1,col))/(2*dZeta);
        dE2v_dE(row,col) = (E2v(row+2,col+1)-E2v(row,col+1))/(2*dEta);
    end
end

% Viscous Zeta Term
for row = 1:1:N
    for col = 1:1:N-1
        zeta_x_hfnd  = (1/2)*(zeta_x(row,col+1)+zeta_x(row,col));
        g11_hfnd     = zeta_x_hfnd^2;
        dv_dZ_hfnd   = (vv(row,col+1)-vv(row,col))/(dZeta);
        J_x_hfnd     = (1/2)*(J(row,col+1)+J(row,col));
        
        E1vv_hfnd(row,col) = (g11_hfnd/(Re*J_x_hfnd))*...
                                (dv_dZ_hfnd);
    end
end

% Viscous Eta Term
for col = 1:1:N
    for row = 1:1:N-1
        eta_y_hfnd   = (1/2)*(eta_y(row+1,col)+eta_y(row,col));
        g22_hfnd     = eta_y_hfnd^2;
        dv_dE_hfnd   = (vv(row+1,col)-vv(row,col))/(dEta);  
        J_y_hfnd     = (1/2)*(J(row+1,col)+J(row,col));
    
        E2vv_hfnd(row,col) = (g22_hfnd/(Re*J_y_hfnd))*...
                                (dv_dE_hfnd);
    end
end

% Dissipation central finite difference
% - D1v is [25x24]
% - D2v is [24x25]
for row = 1:1:N-2
    for col = 1:1:N-2
        dD1v_dZ(row,col) = (D1v(row+1,col+1)-D1v(row+1,col))/dZeta;
    end
end
for row = 1:1:N-2
    for col = 1:1:N-2
        dD2v_dE(row,col) = (D2v(row+1,col+1)-D2v(row,col+1))/dEta;
    end
end
assignin('base','dD1v_dZ',dD1v_dZ);
assignin('base','dD2v_dE',dD2v_dE);

% Viscous Derivatives
% - E1vv_hfnd is [25x24] matrix
% - E2vv_hfnd is [24x25] matrix
for row = 1:1:N-2
    for col = 1:1:N-2
        dE1vv_dZ(row,col) = (E1vv_hfnd(row+1,col+1)-E1vv_hfnd(row+1,col))/dZeta;
        dE2vv_dE(row,col) = (E2vv_hfnd(row+1,col+1)-E2vv_hfnd(row,col+1))/dEta;
    end
end

for row = 1:1:N-2
    for col = 1:1:N-2
        Rv_n(row,col) = J(row+1,col+1)*(-dE1v_dZ(row,col)-dE2v_dE(row,col)+...
                            dE1vv_dZ(row,col)+dE2vv_dE(row,col)-...
                            dD1v_dZ(row,col)-dD2v_dE(row,col));
    end
end

RHSv = Rv_n;

assignin('base','E1v',E1v);
assignin('base','E2v',E2v);
assignin('base','dE1v_dZ',dE1v_dZ);
assignin('base','dE2v_dE',dE2v_dE);



