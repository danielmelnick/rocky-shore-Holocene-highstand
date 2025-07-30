clear
addpath(genpath(pwd))
load HRSL_data.mat

nim=numel(RSL); colorMap = parula(nim);

% Plot
figure(2), clf, %axis off, hold on
h=0.35; w=0.35; bo=.1;

% A: Elevation pdfs
ax1 = axes('position',[1.2*bo 2*bo+h w h]); hold on

for i=1:nim
    plot(RSL(i).pdfeu(:,2),RSL(i).pdfeu(:,1),'color',colorMap(i,:))
end
xlabel('p'), ylabel('Elevation (m)')
ylim([0 7]), xlim([0.02 .55]), 
box on, set(gca,'tickdir','out') %,'XAxisLocation','top')
title('Post-2010 earthquake datum')
% add colorbar
hcb = colorbar('Location','manual','Orientation', 'horizontal','position',[0.3 0.25+h 0.12 0.02]); colormap("parula")
hcb.YTick=[0,1]; hcb.Limits=[0,1]; hcb.LimitsMode='manual';
hcb.YTickLabel={'4','8.4'}; hcb.Label.String={'Age (ka)'};

% B: age-elevation relation
ax2 = axes('position',[1.5*bo+w 2*bo+h w h]); hold on
p3=errorbar([RSL.age],[RSL.pr50],([RSL.pr25]-[RSL.pr75])./2,'.k','markersize',1,'capsize',0,'HandleVisibility','off');
p1=plot([RSL.age],[RSL.pdfmax],'or','MarkerFaceColor','auto');
p2=plot([RSL.age],[RSL.pr50],'ok','MarkerFaceColor','auto');
box on, ylabel('Elevation (m)'), xlabel('Age (ka)')
legend([p1,p2],{'Maximum likelihood','Median'},'location','northeast')
xlim([age1,age2]), ylim([1.5 4.2])
set(gca,'yaxislocation','right')
% add linear regression
p4=plot([RSL.age],f,'--r','HandleVisibility','off'); %,[RSL.age],f+2.*delta,'--r',[RSL.age],f-2.*delta,'--r')
p5=plot([RSL.age],ff,'--k','HandleVisibility','off');

% C: Pre-Maule
ax3 = axes('position',[1.2*bo bo w h]); hold on
for i=1:nim
    plot(RSL2(i).pdfeu(:,2),RSL2(i).pdfeu(:,1),'color',colorMap(i,:))
end
xlabel('p'), ylabel('Elevation (m)')
ylim([0 7]), xlim([0.02 .55]), 
box on, set(gca,'tickdir','out')
title('Pre-2010 earthquake datum')
% add colorbar
hcb = colorbar('Location','manual','Orientation', 'horizontal','position',[0.3 0.15 0.12 0.02]); colormap("parula")
hcb.YTick=[0,1]; hcb.Limits=[0,1]; hcb.LimitsMode='manual';
hcb.YTickLabel={'4','8.4'}; hcb.Label.String={'Age (ka)'}; %axis off

% D: age-elevation relation
ax4 = axes('position',[1.5*bo+w bo w h]); hold on
p3=errorbar([RSL2.age],[RSL2.pr50],([RSL2.pr25]-[RSL2.pr75])./2,'.k','markersize',1,'capsize',0,'HandleVisibility','off');
p1=plot([RSL2.age],[RSL2.pdfmax],'or','MarkerFaceColor','auto');
p2=plot([RSL2.age],[RSL2.pr50],'ok','MarkerFaceColor','auto');
box on, ylabel('Elevation (m)'), xlabel('Age (ka)')
legend([p1,p2],{'Maximum likelihood','Median'},'location','northeast')
xlim([age1,age2]), ylim([1.5 4.2])
set(gca,'yaxislocation','right')

% add linear regression
p4=plot([RSL2.age],f2,'--r','HandleVisibility','off'); %,[RSL.age],f+2.*delta,'--r',[RSL.age],f-2.*delta,'--r')
p5=plot([RSL2.age],ff2,'--k','HandleVisibility','off');

bo=.08; ho=.94;
tx1 = axes('position',[bo ho .1 .1]); axis off
text(0,0,'a','FontSize',12)
tx2 = axes('position',[6*bo ho .1 .1]); axis off
text(0,0,'b','FontSize',12)
tx3 = axes('position',[bo ho/2 .1 .1]); axis off
text(0,0,'c','FontSize',12)
tx4 = axes('position',[6*bo ho/2 .1 .1]); axis off
text(0,0,'d','FontSize',12)

% export
rect = [2,4,20,20]; set(gcf,'PaperUnits','centimeters','PaperType','A4','paperposition',rect,'paperorientation','portrait');
fout = 'figs/HRSL_F3_HHSpdfs.png'; saveas(gcf,fout,'png');
fout = 'figs/HRSL_F3_HHSpdfs.pdf'; %saveas(gcf,fout,'pdf');
opts.Resolution = 150;
opts.BackgroundColor = 'none';
opts.ContentType='vector';
exportPlotToPDF_Advanced(gcf, fout, rect, opts);