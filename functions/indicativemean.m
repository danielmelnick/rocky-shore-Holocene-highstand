function IM = indicativemean(CPD,Lon,Lat,slope)

% Indicative meaning calculator for marine terraces after 
% Lorscheid and Rovere (2019) Open Geospatial Data, Software and Standards 4:10

% INPUT:
% CPD - Coastal Points Database from Lorscheid and Rovere
% Lon - longitude
% Lat - latitude
% slope - beach slope (of empty, average slope of 0.08 will be used)

% OUTPUT:
% IM  - indicative meaning structure with
% IM.MLLW - mean lower low water [m]
% IM.MHHW - mean higher high water [m]
% IM.DB - breaking depth of waves [m]
% IM.Ob - ordinary berm [m]
% IM.Ob2 - max. ordinary berm (storm wave swash height) [m]
% IM.IR - indicative range [m]
% IM.terraceRWL - relative water level of terrace [m]
% IM.platformRWL - relative water level of platform [m]

% nearest point (of CPD) to the input point
k = dsearchn(CPD(:,[7,8]),[Lon Lat]);

% extract values from CPD
Hs = CPD(k,1);                % mean significant wave height [m]
Tp = CPD(k,2);                % mean wave period [s]
Hss = CPD(k,1) + 2*CPD(k,3);  % Hs + 2 standart deviations [m]
Tps = CPD(k,2) + 2*CPD(k,4);  % Tp + 2 standart deviations [s]
IM.MLLW = CPD(k,5);           % mean lower low water [m]
IM.MHHW = CPD(k,6);           % mean higher high water [m]

% beach slope
if isempty(slope) == 1
    S = 0.08;
else
    S = slope;
end

g = 9.81; % [m/s^2]

% deepwater wave length
L0 = g.*(Tp.^2)/(2.*pi);

% breaking depth of waves (lower limit)
IM.DB = (L0*((3.86.*S.^2-1.98.*S+0.88).*((Hs./L0).^0.84)))*(-1);

% storm wave swash height (upper limit)
L02 = g.*(Tps.^2)./(2.*pi);
IM.Ob = 1.1*(0.35*S*(Hs*L0)^0.5 + (0.5*(Hs*L0*(0.563.*S^2+0.004)).^0.5))+IM.MHHW;
IM.Ob2 = 1.1*(0.35*S*(Hss*L02)^0.5 + (0.5*(Hss*L02*(0.563.*S^2+0.004)).^0.5))+IM.MHHW;

% indicative meaning
IM.IR = IM.Ob2 - IM.DB;
IM.terraceRWL = (IM.Ob2 + IM.DB)/2;
IM.platformRWL = (IM.MHHW + (IM.DB + IM.MLLW)/2)/2;
IM.Lat=Lat;
IM.Lon=Lon;

end














