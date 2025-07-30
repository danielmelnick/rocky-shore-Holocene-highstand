clear
addpath(genpath(pwd))
load HRSL_data.mat
load HRSL_giamodels.mat

% get range from min gia to max gia (6.2 to 7.4) post-2010 datum
a=[];
a(:,1)=[RSL.age];
a(:,2)=[RSL.pdfmax];
a(a(:,1)<min([giam.t5]),:)=[];
a(a(:,1)>max([giam.t6]),:)=[];
HRSL(1)=mean(a(:,2));
HRSL(2)=mean(a(:,2))-min(a(:,2));
HRSL(3)=mean(a(:,2))-max(a(:,2));
HRSLp.xi=2:0.01:3.6;
HRSLp.n=normpdf(HRSLp.xi,mean(a(:,2)),std(a(:,2)));

% get range from min gia to max gia (6.2 to 7.4) pre-2010 datum
a=[];
a(:,1)=[RSL2.age];
a(:,2)=[RSL2.pdfmax];
a(a(:,1)<min([giam.t5]),:)=[];
a(a(:,1)>max([giam.t6]),:)=[];
HRSL2(1)=mean(a(:,2));
HRSL2(2)=mean(a(:,2))-min(a(:,2));
HRSL2(3)=mean(a(:,2))-max(a(:,2));
HRSL2p.xi=HRSLp.xi;
HRSL2p.n=normpdf(HRSLp.xi,mean(a(:,2)),std(a(:,2)));

save('plotdata/HRSL.mat',"HRSL","HRSL2")

% age-elevation linear regression at 10yr step
[p,s]=polyfit([RSL.age],[RSL.pdfmax],1);
ti=min([giam.t5]):0.01:max([giam.t6]);
[f,delta] = polyval(p,[RSL.age],s);
[f1,delta1] = polyval(p,ti,s);
[pp,ss]=polyfit([RSL.age],[RSL.pr50],1); 
[ff,ddelta] = polyval(pp,[RSL.age],ss);
[ff1,ddelta1] = polyval(pp,ti,ss);

% calcupate overlap
overlap_area = trapz(HRSLp.xi, min(HRSLp.n, HRSL2p.n));
overlap_percent = overlap_area * 100;

%% get best-fit models and plot 
figure(1), clf, axis off, hold on
h=0.78; 

% plot GIA curves
ax1 = axes('position',[0.1 0.15 0.5 h]); hold on
for k=2:13
    h1=plot([gia.D5i(:,1)],[gia.D5i(:,k)],'-b');
    h2=plot([gia.D6i(:,1)],[gia.D6i(:,k)],'-r');    
end
plot([giam.t5],[giam.m5],'ok',[giam.t6],[giam.m6],'dk','markerfacecolor','k')
h3=plot([RSL.age],f,'-g','linewidth',1.5); %,[RSL.age],f+2.*delta,'--r',[RSL.age],f-2.*delta,'--r')
h4=plot(ti,f1,'-k','linewidth',2.5); %,[RSL.age],f+2.*delta,'--r',[RSL.age],f-2.*delta,'--r')
plot([RSL.age],f+2.*delta,'--g',[RSL.age],f-2.*delta,'--g','linewidth',0.8);
xlabel('Age (ka)'), ylabel('RSL elevation (m)'), box on, set(gca,'tickdir','out')
xlim([4 8.4]), ylim([0 4.6])
Leg=legend([h1,h2,h3,h4],...
    {'ICE-5G','ICE-6G','Post-Maule ML model','Geomorphic estimate'},'location','northwest');

% histograms inset
ax2 = axes('position',[0.15 0.22 0.15 0.2]); hold on; bk=0.9;
set(ax2, 'Color', [bk bk bk])
histogram([giam.t5])
histogram([giam.t6])
xlabel('Age (ka)'), ylabel('n'), box on
xlim([6 7.5]), set(gca,'XTick',6:0.5:7)
ylim([0 7.5])

ax3 = axes('position',[0.31 0.22 0.15 0.2]); hold on
set(ax3, 'Color', [bk bk bk])
histogram([giam.m5])
histogram([giam.m6])
xlabel('Elevation (m)'), ylabel('n'), box on
set(gca,'yaxislocation','right')
xlim([3.3 5.1])
ylim([0 7.5])

% B: Highstand elevations - Litospheric thickness
ax4 = axes('position',[0.65 0.15 0.28 h]); hold on
col = flipud(winter(4)); colormap(ax4,col)
sc=80; 
%S.Vertices = [65 HRSL(1)+HRSL(2); 126 HRSL(1)+HRSL(2); 126 HRSL(1)+HRSL(3); 65 HRSL(1)+HRSL(3)]; S.Faces = [1 2 3 4];
%S.FaceColor = [.8 .8 .8]; %'red';
%S.FaceColor = [0.99, 0.85, 0.13]; S.EdgeColor = [1,1,1]; S.LineWidth = 0.01;
%patch(S)

plot(65:126,HRSL(1)*ones(numel(65:126),1),'-k','linewidth',1.5)
%plot(65:126,HRSL(1)+HRSL(2)*ones(numel(65:126),1),'--k','linewidth',1)
%plot(65:126,HRSL(1)-HRSL(2)*ones(numel(65:126),1),'--k','linewidth',1)
bk=[.7,.7,.7];
plot(65:126,HRSL2(1)*ones(numel(65:126),1),'-','linewidth',1.5,'color',bk)
%plot(65:126,HRSL2(1)+HRSL2(2)*ones(numel(65:126),1),'--b','linewidth',1,'color',bk)
%plot(65:126,HRSL2(1)-HRSL2(2)*ones(numel(65:126),1),'--b','linewidth',1,'color',bk)

h1=scatter([giam.lito]-1,[giam.m5],sc,[giam.visc],'filled','s','MarkerEdgeColor','k');
h2=scatter([giam.lito]+2,[giam.m6],sc,[giam.visc],'filled','d','MarkerEdgeColor','k');
xlabel('Lithosphere thickness (km)'), ylabel('Highstand peak elevation (m)')
set(gca,'yaxislocation','right','YTick',3:0.5:4.5)
legs={'ICE-5G','ICE-6G'}; 
lgd=legend([h1 h2],legs,'location','north'); 
title(lgd,'Ice model')
box on, set(gca,'tickdir','out','xtick',[71,96,120]) %title('ICE 5G')
xlim([66 125]), 
ylim([2.8 4.6])

% Plot PDFs
ax41 = axes('position',[0.65 0.15 0.28 h]); hold on; axis off
fc=0.4*ones(1,3);
fill(HRSLp.n, HRSLp.xi,fc, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
text(2,HRSL(1)+0.06,{'Post-Maule earthquake' 'geomorphic estimate'},'fontsize',9,'Color','k')
fc=0.7*ones(1,3);
fill(HRSL2p.n, HRSL2p.xi,fc , 'FaceAlpha', 0.2, 'EdgeColor', 'none');
text(2,HRSL2(1)-0.05,{'Pre-Maule earthquake' 'geomorphic estimate'},'fontsize',9,'Color',0.7*fc)
ylim([2.8 4.6]), xlim([0.1 4])

ax5 = axes('position',[0.3 1 0.1 0.1]); hold on; axis off
colormap(ax5,col); 
hcb=colorbar('position',[0.66 0.74 0.02 0.18]);
hcb.YTick=[0,0.33,.66,1];
hcb.YTickLabel={'5x10^1^9','8x10^1^9','1x10^2^0','2x10^2^0'};
hcb.Label.String={'Mantle viscosity (Pa*s)'};

bo=.06; ho=.97;
tx1 = axes('position',[bo ho .1 .1]); axis off
text(0,0,'a','FontSize',12)
tx2 = axes('position',[10.5*bo ho .1 .1]); axis off
text(0,0,'b','FontSize',12)

%% export
rect = [2,4,25,16]; 
set(gcf,'PaperUnits','centimeters','PaperType','A4','paperposition',rect,'paperorientation','portrait');
fout = 'figs/HRSL_F4_GIA_HRSL.png'; saveas(gcf,fout,'png');
fout = 'figs/HRSL_F4_GIA_HRSL.pdf'; %saveas(gcf,fout,'pdf');
opts.Resolution = 150;
opts.BackgroundColor = 'none';
opts.ContentType='vector';
exportPlotToPDF_Advanced(gcf, fout, rect, opts);
