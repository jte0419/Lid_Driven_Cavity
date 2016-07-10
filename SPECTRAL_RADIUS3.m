function [pA1,pA2] = SPECTRAL_RADIUS3(U,V,g11,g22,J,N)

% INPUTS
% - U : Contravariant velocity in zeta direction [NxN]
% - V : Contravariant velocity in eta direction [NxN]
% - g : Metric tensor of the transformation [Nx1]
% - J : Jacobian of the transformation [NxN]
%
% OUTPUTS
% - specRad1 : Spectral radius matrix in zeta direction
% - specRad2 : Spectral radius matrix in eta direction

% Compute spectral radius pA1
% - X or Zeta direction
oneU = abs(U);
twoU = sqrt((U.^2)+g11);
pA1 = 1./((oneU+twoU).*J);

% Calculate spectral radius pA2
% - Y or Eta direction
oneV = abs(V);
twoV = sqrt((V.^2)+g22);
pA2 = 1./((oneV+twoV).*J);