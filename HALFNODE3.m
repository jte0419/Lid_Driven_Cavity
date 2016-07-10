function [pA1_hfnd,pA2_hfnd] = HALFNODE3(pA1,pA2,N)

% INPUTS
% - pA1 : Spectral radius in the zeta direction [NxN]
% - pA2 : Spectral radius in the eta direction [NxN]
% - N   : Number of grid nodes
%
% OUTPUTS
% - pA1_hfnd : Spectral radius at the half nodes [N-1xN-1]
% - pA2_hfnd : Spectral radius at the half nodes [N-1xN-1]

% Pre-allocate matrices
pA1_hfnd = zeros(N,N-1);    % One less node in x direction
pA2_hfnd = zeros(N-1,N);    % One less node in y direction

% Halfnode values - Zeta
% - [25x24] matrix
for row = 1:1:N
    for col = 1:1:N-1
        % Differencing left to right in matrix
        pA1_hfnd(row,col) = (pA1(row,col+1)-pA1(row,col));
    end
end

% Halfnode values - Eta
% - [24x25] matrix
for row = 1:1:N-1
    for col = 1:1:N
        % Differencing top to bottom in matrix
        pA2_hfnd(row,col) = (pA2(row+1,col)-pA2(row,col));
    end
end

