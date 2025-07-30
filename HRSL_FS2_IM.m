% plot Figure S3 - indicative meaning histogram

clear
addpath(genpath(pwd))
load IM_HRSL.mat

figure(1), clf

h=0.7; w=0.7;
ax1 = axes('position',[0.2 0.2 w h]); hold on

histogram([IM.terraceRWL]);
ylabel('Number of sites'), xlabel('Elevation (m)')
title('Marine terrace indicative meaning (34{\circ}-38{\circ}S)')
box on, set(gca,'tickdir','out')
text(-0.15,60,sprintf('n=%u sites\nMean=%1.2f ± %1.2f m', ...
    numel(IM),mean([IM.terraceRWL]),std([IM.terraceRWL])),'FontSize',10)

% export
rect = [2,4,15,15]; 
set(gcf,'PaperUnits','centimeters','PaperType','A4','paperposition',rect,'paperorientation','portrait');
fout = 'figs/HRSL_FS2_IM.png'; saveas(gcf,fout,'png');
%fout = 'figs/HRSL_FS3_IM.pdf'; %saveas(gcf,fout,'pdf');
%opts.Resolution = 150; opts.BackgroundColor = 'none'; opts.ContentType='vector';
%exportPlotToPDF_Advanced(gcf, fout, rect, opts);
