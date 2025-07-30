%% plot Figure S2 - GIA model curves
clear
addpath(genpath(pwd))
load HRSL_giamodels.mat

D5i=gia.D5i;
D6i=gia.D6i;

figure(1), clf

% 5G
h=0.6; w=0.4;
ax1 = axes('position',[0.1 0.2 w h]); hold on
for j=1:4        
     h1=plot(D5i(:,1),D5i(:,1+j),'-','color',[0 j/4 1-j/4]);
     h2=plot(D5i(:,1),D5i(:,1+4+j),'--','color',[0 j/4 1-j/4]);
     h3=plot(D5i(:,1),D5i(:,1+8+j),'-.','color',[0 j/4 1-j/4]);
end
ylim([0 4.5]), xlim([0 9])
box on, set(gca,'tickdir','out'), title('ICE 5G')
xlabel('Time (ka)'), ylabel('Elevation (m)')
plot([giam.t5],[giam.m5],'ok')

% colormap
for j=1:64
    col(j,1)=0; col(j,2)=j/64; col(j,3)=1-j/64;
end
col=flipud(col);
% colorbar
ax2 = axes('position',[0.2 1 0.1 0.1]); hold on
colormap(col); 
hcb=colorbar('position',[0.2 0.5 0.02 0.25]);
hcb.YTick=[0,0.33,.66,1];
hcb.YTickLabel={'5x10^1^9','8x10^1^9','1x10^2^0','2x10^2^0'};
hcb.Label.String={'Mantle viscosity (Pa*s)'};
axis off

% 6G
ax1 = axes('position',[0.12+w 0.2 w h]); hold on
for j=1:4        
     h1=plot(D6i(:,1),D6i(:,1+j),'-','color',[0 j/4 1-j/4]);
     h2=plot(D6i(:,1),D6i(:,1+4+j),'--','color',[0 j/4 1-j/4]);
     h3=plot(D6i(:,1),D6i(:,1+8+j),'-.','color',[0 j/4 1-j/4]);
end
ylim([0 4.5]), xlim([0 9])
box on, set(gca,'tickdir','out','yaxislocation','right'), title('ICE 6G')
xlabel('Time (ka)'), ylabel('Elevation (m)')
plot([giam.t6],[giam.m6],'dk')

legs={'71 km','96 km','120 km'}; lgd=legend(legs,'location','northwest'); title(lgd,'Lithosphere thickness')
xlabel('Time (ka)'), ylabel('Elevation (m)')

% export
rect = [2,4,18,14]; 
set(gcf,'PaperUnits','centimeters','PaperType','A4','paperposition',rect,'paperorientation','portrait');
fout = 'figs/HRSL_FS1_GIAmodels.png'; saveas(gcf,fout,'png');
%fout = 'figs/HRSL_FS1_GIAmodels.pdf'; %saveas(gcf,fout,'pdf');
%opts.Resolution = 150; opts.BackgroundColor = 'none'; opts.ContentType='vector';
%exportPlotToPDF_Advanced(gcf, fout, rect, opts);
