function [PBC,uBC,vBC] = BOUNDARY_CONDITIONS3(PBC,uBC,vBC,velLid,N)

% INPUTS
% - PBC    : Pressure [N-2xN-2]
% - uBC    : Velocity in x/zeta direction [N-2xN-2]
% - vBC    : Velocity in y/eta direction [N-2xN-2]
% - velLid : Velocity of the lid [m/s]
% - N      : Number of grid nodes
% 
% OUTPUTS
% - PBC: Pressure with boundary conditions [NxN]
% - uBC: X/Zeta velocity with boundary conditions [NxN]
% - vBC: Y/Eta velocity with boundary conditions [NxN]

% U Velocity Boundary Conditions
uBC(:,N-1) = 0;                     % Right wall
uBC(:,N)   = 0;                     % Left wall
uBC(N-1,:) = 0;                     % Bottom wall
uBC(N,:)   = velLid;                % Top wall (lid)
uBC = circshift(uBC,[1 1]);         % Shift matrix once in x and y direction

% V Velocity Boundary Conditions
vBC(N-1,:) = 0;                     % Bottom wall
vBC(N,:)   = 0;                     % Top wall (lid)
vBC(:,N-1) = 0;                     % Right wall
vBC(:,N)   = 0;                     % Left wall
vBC = circshift(vBC,[1 1]);         % Shift matrix once in x and y direction

% Pressure Boundary Conditions
PBC(N-1,:) = PBC(N-2,:);            % Bottom wall
PBC(N,:)   = PBC(1,:);              % Top wall (lid)
PBC(:,N-1) = PBC(:,N-2);            % Right wall
PBC(:,N)   = PBC(:,1);              % Left wall
PBC = circshift(PBC,[1 1]);         % Shift matrix once in x and y direction


