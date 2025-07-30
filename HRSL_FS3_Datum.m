clear
addpath(genpath(pwd))

% load Terrace data
S=shaperead('HRSL_terracedata.shp');
% Remove MIS5 sites
S([S.level]<5)=[]; 

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
    S(i).Zpre=S(i).Z-S(i).U2010;    
end

% Export data
shapewrite(S,'indata/HRSL_terracedata_pre2010.shp')

% Plot
figure(1), clf, axis off
ax = axes('position',[0.1 0.1 0.8 0.8]); hold on; bk=0.9;
set(ax, 'Color', [bk bk bk])
plot(2:14,2:14,'-k'), box on
scatter([S.Z],[S.Zpre],50,[S.lat],'filled','MarkerEdgeColor','w')
xlabel('Post-2010 earthquake elevation (m)'), ylabel('pre-2010 earthquake elevation (m)')
hcb=colorbar('position',[0.23 0.5 0.03 0.3]);
hcb.Label.String={'Latitude along coast'};
hcb.TickDirection="out";
% Export fig
rect=[0,0,15,15]; %[xmin ymin width height]
set(gcf,'PaperType','A4','PaperUnits','centimeters','Paperposition',rect);
fout='figs/HRSL_FS3_corr_27F.png'; saveas(gca,fout,'png')
%%
fout='figs/HRSL_FS4_corr_27F.pdf'; %saveas(gca,fout,'png')
opts.Resolution = 150;
opts.BackgroundColor = 'none';
opts.ContentType='vector';
exportPlotToPDF_Advanced(gcf, fout, rect, opts);






