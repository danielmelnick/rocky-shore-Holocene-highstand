clear
addpath(genpath(pwd))
load HRSL_data.mat
load MIS1.mat
load GIAM5_HRSL_LEM_postM.mat
load GIAM6_HRSL_LEM_postM.mat
load HRSL_LEM_Best_postM.mat

%% plot GIA model misfits
figure(1), clf, h=0.6; w=0.2; sz=10; fz=11;

ax1 = axes('position',[0.1 0.2 w h]); hold on
j=12;
hh=plot([GIAM5(j).L]-1,[GIAM5(j).dshz],'ok','MarkerSize',25);          
set(get(get(hh,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
col = flipud(winter(4));
for j=1:12     
    k=[GIAM5(j).V];
    h1=plot([GIAM5(j).L]-1,[GIAM5(j).dshz],'s','markerfacecolor',col(k,:),'MarkerSize',sz,'MarkerEdgeColor','k');          
    h2=plot([GIAM6(j).L]+1,[GIAM6(j).dshz],'d','markerfacecolor',col(k,:),'MarkerSize',sz,'MarkerEdgeColor','k');          
end
xlim([65 126])
box on, set(gca,'tickdir','out','xtick',[71;96;120])
xlabel('Lithosphere thickness (km)','fontsize',fz), ylabel('Mean LEM - TerraceM shoreline elevation (m)','fontsize',fz)
legs={'ICE-5G','ICE-6G'}; lgd=legend(legs,'location','northwest'); %title(lgd,'Lithosphere thickness')

% colorbar
ax2 = axes('position',[0.05 1 0.1 0.1]); hold on
colormap(ax2,col);
hcb=colorbar('position',[0.35 0.55 0.025 0.2]);
hcb.YTick=[0,0.33,.66,1];
hcb.YTickLabel={'5x10^1^9','8x10^1^9','1x10^2^0','2x10^2^0'};
hcb.Label.String={'Mantle viscosity (Pa*s)'};
axis off

MHH=0; h=0.6; w=0.6; sz=10; xma=16;

% LEM-TerraceM shoreline elevation
k=1; LEMB=mod(k).LEM;

ax3 = axes('position',[0.35 0.2 w h]); hold on; bk=0.9;
set(ax3, 'Color', [bk bk bk])

col2 = parula(numel(LEMB));
colormap(ax3,col2)
plot(0:xma,0:xma,'-k')
errorbar([LEMB.shz]-MHH,[LEMB.MIS1z],[mis1.Ze],'dk','CapSize',0)
scatter([LEMB.shz]-MHH,[LEMB.MIS1z],50,[LEMB.UR],'d','filled','MarkerEdgeColor','w')
ylabel('TerraceM shoreline elevation (m)','fontsize',fz), xlabel('LEM shoreline elevation (m)','fontsize',fz), box on
title('Best-fit GIA sea-level curve')
R=corrcoef([LEMB.shz],[LEMB.MIS1z]);
text(6,14,sprintf('R^2=%2.2f\nn=%u',R(1,2)^2,numel(mis1)),'FontSize',12)
axis([2 xma 2 xma]), box on, 
axis equal square
set(gca,'YAxisLocation','right','TickDir','out')

% colorbar
ax5 = axes('position',[0.05 1 0.1 0.1]); hold on; axis off
colormap(ax5,col2);
hcb=colorbar('position',[0.35 0.2 0.025 0.2]);
hcb.Label.String={'Uplift rate (mm/yr)'}; 
hcb.LimitsMode='manual';
hcb.Limits=[min([LEMB.UR]),max([LEMB.UR])];
hcb.YTick=0.5:0.5:1.5;

ax4 = axes('position',[0.7 0.28 .1 .13]); hold on
histogram([LEMB.dshz]) 
set(gca,'YAxisLocation','right')
xlim([-2 3.5])
ylabel('n'), xlabel('Residual elevation (m)'), box on
title(sprintf('Mean=%3.2f m',mean([LEMB.dshz])),'FontSize',10)
bo=.05; ho=.85;
tx1 = axes('position',[bo ho .1 .1]); axis off
text(0,0,'a','FontSize',12)
tx2 = axes('position',[bo+w/1.5 ho .1 .1]); axis off
text(0,0,'b','FontSize',12)

rect = [2,4,25,15]; 
set(gcf,'PaperUnits','centimeters','PaperType','A4','paperposition',rect,'paperorientation','portrait');
fout = 'figs/HRSL_F5_LEM.png'; saveas(gcf,fout,'png');
fout = 'figs/HRSL_F5_LEM.pdf'; %saveas(gcf,fout,'pdf');
opts.Resolution = 150;
opts.BackgroundColor = 'none';
opts.ContentType='vector';
exportPlotToPDF_Advanced(gcf, fout, rect, opts);

