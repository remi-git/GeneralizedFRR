function [TEA_Cxx] = Cross_TA_Cxx_flip(X, Y, MaxT_0, NTime, NTraj)
% FLIP for inverse time and ensemle ordering
% Computes time average Cxx
%   LagArray is the time in SECONDS of each individual trajectory
EaCxx = zeros(NTraj, MaxT_0);
Wb = waitbar(0, 'correlation');
for i = 1:NTraj % mixed time-ensemble averaging
    waitbar(i/(NTraj))
    Cxx_in = zeros(NTime-MaxT_0, MaxT_0);
    for t_0 = 1:NTime-MaxT_0
        Cxx_in(t_0, :) = X(t_0, i) .* Y(t_0:t_0+MaxT_0-1, i);
    end
    EaCxx(i,:) = mean(Cxx_in, 1);
end
TEA_Cxx = mean(EaCxx, 1);
delete(Wb)
end