clear
addpath(genpath(pwd))

% load Terrace data
S=shaperead('HRSL_terracedata.shp');

% Load and grid coseismic
U=load('Moreno2012_coseismic_uplift.txt');
U(U(:,2)<-39,:)=[];
U(U(:,1)>-72,:)=[];
F = scatteredInterpolant(U(:,1),U(:,2),U(:,3));
F.Method = 'natural';

% interpolate coseismic to terrace sites
for i=1:numel(S)
    [la,lo]=utm18_2deg(S(i).X,S(i).Y);
    S(i).lat=la;
    S(i).lon=lo;
    S(i).U2010=F(lo,la);
    if S(i).level==5        
        S(i).Zpre=S(i).Z-S(i).U2010;    
    else
        S(i).Zpre=S(i).Z;
    end
end
% Export data
shapewrite(S,'indata/HRSL_terracedata.shp')







