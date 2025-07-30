function S = S_upliftrate2(S,Z,SL)
%function [H,U,Ue] = S_upliftrate2(S,SL)

% Calculate uplift rate from shoreline angle elevation
% using the formula from Gallen et al. (2018) EPSL

% inputs
% E:  terrace elevation
% dE: terrace elevation error
% e:  RSL elevation
% de: RSL elevation error
% T:  RSL time
% dT: RSL time error
% Z:  Name of the elevation field to use

%S.H=[S.Zrsl]-[SL.e];
S.H=[Z]-[SL.e];
dH=S.Ze^2+SL.de^2;

S.U = [S.H]/[SL.T];
S.Ue=[S.U]^2*((dH^2/[S.H]^2) + ([SL.dT]^2/[SL.T]^2));
