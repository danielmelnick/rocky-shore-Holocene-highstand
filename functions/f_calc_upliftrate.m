function [U,Ue] = f_calc_upliftrate(Z,Ze,SL)

% Calculate uplift rate from shoreline angle elevation
% using the formula from Gallen et al. (2018) EPSL

% inputs
% Z:  terrace elevation
% Ze: terrace elevation error
% SL.e:  RSL elevation
% SL.de: RSL elevation error
% SL.T:  RSL time
% SL.dT: RSL time error

%S.H=[S.Zrsl]-[SL.e];
H=Z-SL.e;
dH=Ze.^2+SL.de.^2;

U = H./SL.T;
Ue=U.^2*((dH.^2/H.^2) + ([SL.dT].^2/[SL.T].^2));

%S.U=U;
%S.Ue=Ue;

end