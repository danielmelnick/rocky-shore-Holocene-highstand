clear
addpath(genpath(pwd))
load HRSL_data.mat
clear line
TE=shaperead("PB2002_boundaries.shp");
load coast.mat

mis5e(92)=[]; %remove one high error site

% Plot
hFig=figure(1); clf, axis off
set(hFig, 'Color', 'w')

bo=0.1; h=0.38; w=0.8; m=4; mks=8; mw=8; m5col=[0.2,0.72,0.63]; fs=8; bk=0.9;

% A - Map
hax1 = axes('position',[bo/3 .51 .32 .32]); hold on
set(hax1, 'Color', [0.75 0.85 0.95])
colormap(hax1,"copper")
latlim = [-39 -33]; lonlim = [-74.5 -71.5];
geoshow(L1, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor','k')
scatter([mis1.Lon],[mis1.Lat],15,[mis1.Zrsl],"filled",'o','MarkerEdgeColor','w')
% Plot Trench
plot([TE.X],[TE.Y],'-r','LineWidth',1.5)
% Add plate convergence vector
lat=-34; lon=-74; vl = plate_vel(lat, lon, -25.4, -124.6, 0.11);    %vl = [vn, ve, vd]
quiver(lon,lat,vl(1,1),-vl(2,1),0.08,'linewidth',2,'color','k','MaxHeadSize',3)
text(-74.2,-33.5,'Nazca'), text(-72.8,-33.6,'SAM')
axis equal, ylim([-39 -33]), xlim([-74.5 -71.5]), box on 
set(gca,'XTick',-74:2:-72,"XTickLabel",{'74^{\circ}W','72^{\circ}W'},...
    'YTick',-38:2:-34,"YTickLabel",{'38^{\circ}S','36^{\circ}S','34^{\circ}S'},'XAxisLocation','top','FontSize',8,'TickDir','out')
xlabel('Latitude'), ylabel('Longitude')
%colorbar
hcbm=colorbar('position',[0.21 .54 0.015 0.07]);
hcbm.Label.String ='Holocene elevation (m)';
hcbm.Label.HorizontalAlignment="center";
pos = hcbm.Label.Position;
pos(2) = pos(2) + 3.9;  
pos(1) = pos(1) + 0.6;  
hcbm.Label.Position = pos;
hcbm.FontSize = 7.2;
hcbm.TickDirection='out';
hcbm.Ticks=[4,8,12];

% B - Along-coast profile 
hax2 = axes('position',[bo bo/1.5 w h]); hold on
set(hax2, 'Color', [bk bk bk])

% Plot MIS-5 elevations
yyaxis left; hold on
plot([mis5e.d_profile]./1e3,[mis5e.U],'.-b','markersize',mks)
he=errorbar([mis5e.d_profile]./1e3,[mis5e.U],[mis5e.Ue],'.b'); he.CapSize=0;
plot([mis5c.d_profile]./1e3,[mis5c.U],'.','markersize',mks,'color',m5col)
he=errorbar([mis5c.d_profile]./1e3,[mis5c.U],[mis5c.Ue],'.','Color',m5col); he.CapSize=0;
xlim([0 500]), set(gca,'tickdir','out'), box off 
ylim([0 2])
ylabel('MIS-5 uplift rate (mm/yr)','fontsize',fs)
set(hax2, 'YColor', 'b');
line([0,500],[2,2],'Marker','none','LineWidth',1,'Color','k')

% Plot Holocene elevations
yyaxis right, hold on
plot([mis1.d_profile]./1e3,[mis1.Zrsl],'.-k','markersize',mks)
he=errorbar([mis1.d_profile]./1e3,[mis1.Zrsl],[mis1.Ze],'.k'); he.CapSize=0;
ylabel('Holocene elevation (m)','fontsize',fs)
xlabel('Distance along profile (km)','fontsize',fs)
xlim([0 500]), set(gca,'tickdir','out'), box off
ylim([0 14])
text(200,9,sprintf('MIS-5c (n=%u)',numel(mis5c)),'fontsize',10,'color',m5col)
text(200,10.5,sprintf('MIS-5e (n=%u)',numel(mis5e)),'color','b','fontsize',10)
text(200,12,sprintf('Holocene (n=%u)',numel(mis1)),'color','k','fontsize',10)
h_ax2 = gca;
h_ax2.YColor = 'k'; 
h_ax2.FontSize = fs;

% C - MIS5-Holocene elevation relation

bo=0.1; h=0.25; w=h; m=4; mks=8; mw=8; xpp=2.3*h;
xbar=0.55; ybar=0.53; scatz=20;

hax3 = axes('position',[1.5*w xpp w h]); hold on
set(hax3, 'Color', [bk bk bk], 'XColor', 'b');
scatter([mis1.mis5Zrsl],[mis1.Z],scatz,[mis1.d_profile]/1e3,'filled','MarkerEdgeColor','w') 
% add regression
d1=[]; d1(:,1)=[mis1.mis5Zrsl]; d1(:,2)=[mis1.Z]; 
d1(isnan(d1(:,2)),:)=[]; d1(isnan(d1(:,1)),:)=[]; d1=sortrows(d1,1);
[p1,s1]=polyfit(d1(:,1),d1(:,2),1); [p1_do,~]=polyval(p1,d1(:,1),s1);
r=corrcoef(d1(:,1),d1(:,2)); r2=r(1,2)^2;
plot(d1(:,1),p1_do,'--r','linewidth',1), 
text(130,4,sprintf('R^2=%3.2g\nn=%u',r2,numel(mis1)),'fontsize',10,'color','k')
box on, ylim([2 13]), xlim([-10 240])
ylabel('Holocene elevation (m)','fontsize',fs);
xlabel('MIS-5e elevation (m)','fontsize',fs);
set(gca,'fontsize',fs,'tickdir','out','xaxislocation','top')

hcb = colorbar('Location','manual','Orientation', 'horizontal','position',[xbar ybar 0.2 0.02]);
hcb.AxisLocationMode = 'manual';
hcb.Label.String ='Distance along profile (km)';
hcb.TickDirection='out';
hcb.Ticks=100:100:400;
hcb.AxisLocation = "in";
%hcb.Label.Position = [mean(hcb.Limits), 1];

% D - Holocene / MIS5e&5c Uplift rate regression
hax4 = axes('position',[2.6*w xpp w h]); hold on, fs=8;
set(hax4, 'Color', [bk bk bk], 'XColor', 'b');

d1=[]; d1(:,1)=[mis5.U]; d1(:,2)=[mis5_1]; 
d1(isnan(d1(:,1)),:)=[]; d1(isnan(d1(:,2)),:)=[]; 
d1=sortrows(d1,1); [p1,s1]=polyfit(d1(:,1),d1(:,2),1); [p1_do,~]=polyval(p1,d1(:,1),s1);
r=corrcoef(d1(:,1),d1(:,2)); r2=r(1,2)^2;
scatter([mis5.U],[mis5_1],scatz,[mis5.Lat],'filled','o','MarkerEdgeColor','w')
plot(d1(:,1),p1_do,'--r','linewidth',1), 
text(1,4,sprintf('R^2=%3.2g\nn=%u',r2,numel(mis5)),'fontsize',10)
box on, ylim([2 13]), xlim([-0.1 1.9])
ylabel('Holocene elevation (m)','fontsize',fs);
xlabel('MIS-5e&5c uplift rate (mm/yr)','fontsize',fs);
set(hax4,'fontsize',fs,'tickdir','out','xaxislocation','top','yaxislocation','right')

% scalebar

%hcb=colorbar('position',[xbar ybar 0.15 0.02]);
%hcb.Orientation = 'horizontal'; hcb.Label.String ='Latitude'; hcb.TickDirection='out'; hcb.Ticks=[-37,-35];
%hcb.TickLabels={'37^{\circ}S','35^{\circ}S'};
%hcb.AxisLocationMode='manual';

% figure panel labels
bo=.04; ho=.9;
tx1 = axes('position',[bo ho .1 .1]); axis off
text(0,0,'a','FontSize',12)
tx2 = axes('position',[bo .5 .1 .1]); axis off
text(0,0,'b','FontSize',12)
tx3 = axes('position',[2*bo+w ho .1 .1]); axis off
text(0,0,'c','FontSize',12)
tx4 = axes('position',[bo+2.35*w ho .1 .1]); axis off
text(0,0,'d','FontSize',12)

% export
rect = [2,4,17,17]; set(gcf,'PaperUnits','centimeters','PaperType','A4','paperposition',rect);
fout = 'figs/HRSL_F1_profile.png'; saveas(gcf,fout,'png');
%%
fout = 'figs/HRSL_F1_profile.pdf'; 
opts.Resolution = 150;
%opts.BackgroundColor = 'none';
opts.ContentType='vector';
exportPlotToPDF_Advanced(gcf, fout, rect, opts);

