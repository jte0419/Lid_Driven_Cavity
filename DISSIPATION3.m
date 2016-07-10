function [D1c,D1u,D1v,D2c,D2u,D2v] = DISSIPATION3(eps,pA1,pA2,P,u,v,N)

% INPUTS
% - eps : Dissipation constant
% - pA1 : Spectral radius in zeta direction [NxN-1]
% - pA2 : Spectral radius in eta direction [N-1xN]
% - P   : Pressure field [NxN]
% - u   : U velocity field [NxN]
% - v   : V velocity field [NxN]
% - N   : Number of grid nodes
% 
% OUTPUTS
% - D1c : Zeta direction dissipation for continuity
% - D1u : Zeta direction dissipation for x-momentum
% - D1v : Zeta direction dissipation for y-momentum
% - D2c : Eta direction dissipation for continuity
% - D2u : Eta direction dissipation for x-momentum
% - D2v : Eta direction dissipation for y-momentum

D1c = zeros(N,N-1);
D1u = zeros(N,N-1);
D1v = zeros(N,N-1);
D2c = zeros(N-1,N);
D2u = zeros(N-1,N);
D2v = zeros(N-1,N);

% Calculate dissipation in zeta (1) direction
for row = 1:1:N
    for col = 2:1:N-2
        D1c(row,col) = eps*pA1(row,col)*(P(row,col+2)-3*P(row,col+1)+3*P(row,col)-P(row,col-1));
        D1u(row,col) = eps*pA1(row,col)*(u(row,col+2)-3*u(row,col+1)+3*u(row,col)-u(row,col-1));
        D1v(row,col) = eps*pA1(row,col)*(v(row,col+2)-3*v(row,col+1)+3*v(row,col)-v(row,col-1));
    end
end

% Set left boundary value to next inner node
D1c(:,1) = D1c(:,2);
D1u(:,1) = D1u(:,2);
D1v(:,1) = D1v(:,2);

% Set right boundary value to next inner node
D1c(:,N-2) = D1c(:,N-1);
D1u(:,N-2) = D1u(:,N-1);
D1v(:,N-2) = D1v(:,N-1);

% Calculate dissipation in eta (2) direction
for row = 2:1:N-2
    for col = 1:1:N
        D2c(row,col) = eps*pA2(row,col)*(P(row+2,col)-3*P(row+1,col)+3*P(row,col)-P(row-1,col));
        D2u(row,col) = eps*pA2(row,col)*(u(row+2,col)-3*u(row+1,col)+3*u(row,col)-u(row-1,col));
        D2v(row,col) = eps*pA2(row,col)*(v(row+2,col)-3*v(row+1,col)+3*v(row,col)-v(row-1,col));
    end
end

% Set top boundary value to next inner node
D2c(1,:) = D2c(2,:);
D2u(1,:) = D2u(2,:);
D2v(1,:) = D2v(2,:);

% Set bottom boundary value to next inner node
D2c(N-1,:) = D2c(N-2,:);
D2u(N-1,:) = D2u(N-2,:);
D2v(N-1,:) = D2v(N-2,:);


