% Script to import Holocene and MIS-5 shorelines, project points along coast-parallel profile, 
% interpolate MIS-5 uplift rates at Holocene locations, and calculate Holocene sea-level elevation 
% D.Melnick 2025

clear
addpath(genpath(pwd))

% Load Holocene TerraceM shorelines
file='indata/HRSL_terracedata.shp'; F=shaperead(file);

% Project points to a trench-paralel profile
line='indata/profile_line_u18S.shp'; 
Pf = shaperead(line); [PLat,~] = utm18_2deg(Pf.X,Pf.Y);
fout='indata/HRSL_terracedata_projected.shp';
S=points2profile(line,file,fout); 

ns=numel(S);
for i=1:ns
   [lat,lon]=utm18_2deg(S(i).X,S(i).Y);
   S(i).Lat=lat;
   S(i).Lon=lon;
   S(i).PLat=S(i).d_profile_Lat+PLat(1);
end

% Sea level data from Rohling et al. 2009 
SL(1).MIS='5e';
SL(1).e=5;  % elevation
SL(1).de=2; % elevation uncertainty
SL(1).T=125;% age
SL(1).dT=5; % age uncertainty
SL(1).ref='Rohling2009'; % source

SL(2).MIS='5c';
SL(2).e=-21;
SL(2).de=2;
SL(2).T=105;
SL(2).dT=5;
SL(2).ref='Rohling2009';

% calculate Indicative Meaning and RSL elevation
load CPD.mat
for i=1:ns
    IM(i) = indicativemean(CPD,[S(i).Lon],[S(i).Lat],0.08);
    S(i).IM = IM(i).terraceRWL;
    S(i).Zrsl = S(i).Z - S(i).IM; %post 2010 datum
    S(i).Zrsl_pre = S(i).Zpre - S(i).IM; % pre 2010 datum
end
save("plotdata/IM_HRSL.mat","IM")

% calculate uplift rates and propagate errors
for i=1:ns
    if S(i).level==1 %5e
        [U,Ue] = f_calc_upliftrate(S(i).Zrsl,S(i).Ze,SL(1));
        S(i).U = U;
        S(i).Ue = Ue;
   elseif S(i).level==4 %5c               
        [U,Ue] = f_calc_upliftrate(S(i).Zrsl,S(i).Ze,SL(2));        
        S(i).U = U;
        S(i).Ue = Ue;
   elseif S(i).level==5 %Holocene        
        S(i).U=[];
        S(i).Ue=[];
   end
end
[S,~] = sortStruct2(S, 'level'); 

% Separate structures for easier handling
% S.shoreline fields: 1=mis5e; 4=mis5c; 5=Holocene
mis5e=S([S.level]==1);
mis5c=S([S.level]==4);
mis5=vertcat([mis5c],[mis5e]);
mis1=S([S.level]==5); 

% Sort along profile
[mis1,~] = sortStruct2(mis1, 'd_profile'); 
[mis5e,~] = sortStruct2(mis5e, 'd_profile'); 
[mis5c,~] = sortStruct2(mis5c, 'd_profile'); 
[mis5,~] = sortStruct2(mis5, 'd_profile'); 

save('plotdata/MIS1.mat','mis1')

% Interpolate MIS5 values to MIS1 sites
mis1a=interp1([mis5e.d_profile],[mis5e.Zrsl],[mis1.d_profile]);
mis1b=interp1([mis5e.d_profile],[mis5e.U],[mis1.d_profile]);
mis1c=interp1([mis5c.d_profile],[mis5c.U],[mis1.d_profile]);

% Interpolate MIS1 elevations (minus IM) to MIS5e and 5c sites
mis5_1=interp1([mis1.d_profile],[mis1.Zrsl],[mis5.d_profile]);
mis5_1_pre=interp1([mis1.d_profile],[mis1.Zrsl_pre],[mis5.d_profile]);

% Add interpolated values to structure
for i=1:numel(mis1)
    mis1(i).mis5Zrsl=mis1a(i);    
    mis1(i).mis5eU=mis1b(i);    
    mis1(i).mis5cU=mis1c(i);        
end

% build coast-parallel interpolated profile
method='Linear';

prof.d=0:0.1:500; %km
prof.Lat=interp1([mis5.d_profile]./1e3,[mis5.Lat],[prof.d],method);
prof.mis5ez=interp1([mis5e.d_profile]./1e3,[mis5e.Z],[prof.d],method);
prof.mis5eu=interp1([mis5e.d_profile]./1e3,[mis5e.U],[prof.d],method);
prof.mis5eue=interp1([mis5e.d_profile]./1e3,[mis5e.Ue],[prof.d],method);
prof.mis5cz=interp1([mis5c.d_profile]./1e3,[mis5c.Z],[prof.d],method);
prof.mis5cu=interp1([mis5c.d_profile]./1e3,[mis5c.U],[prof.d],method);
prof.mis5cue=interp1([mis5c.d_profile]./1e3,[mis5c.Ue],[prof.d],method);
prof.mis1z=interp1([mis1.d_profile]./1e3,[mis1.Z],[prof.d],method);

%% Build RSL history 
i=1; clear RSL, clear B
age1=4; 
age2=8.4; 
dage=0.2;

% use mis5e uplift rate to remove tectonic uplift from mis1 elevations
for age=age1:dage:age2
    
    % tectonic uplift in m for Highstand age with MIS5 rate
    mis1Zpred=[mis1.mis5eU] * age;

    % measured elev. - tectonic uplift % Post2010 datum
    dU=[mis1.Zrsl] - mis1Zpred;
    
    RSL(i,1).age=age;
    RSL(i,1).eustatic=dU;
    RSL(i,1).meaneu=nanmean(dU);     
    pr = prctile(dU,[25 50 75]);
    RSL(i,1).pr25=pr(1);
    RSL(i,1).pr50=pr(2);
    RSL(i,1).pr75=pr(3);    
    [f,xi] = ksdensity(dU);     
    [~, index_to_max] = max(f);
    
    %finds max likelyhood
    RSL(i,1).pdfmax(:,1)=xi(index_to_max);
    RSL(i,1).pdfmaxp(:,1)=f(index_to_max);    
    RSL(i,1).pdfeu(:,1)=xi;
    RSL(i,1).pdfeu(:,2)=f;
    %[phat,pci] = mle([xi,f]);
    %RSL(i,1).MLE_phat=phat;
    %RSL(i,1).MLE_pci=pci;
    %B(:,i)=dU;
    
    % measured elev. - tectonic uplift %Pre-2010 datum
    dU2=[mis1.Zrsl_pre] - mis1Zpred;
    
    RSL2(i,1).age=age;
    RSL2(i,1).eustatic=dU2;
    RSL2(i,1).meaneu=nanmean(dU2);     

    pr = prctile(dU2,[25 50 75]);
    RSL2(i,1).pr25=pr(1);
    RSL2(i,1).pr50=pr(2);
    RSL2(i,1).pr75=pr(3);    

    [f,xi] = ksdensity(dU2);     
    [max_f, index_to_max] = max(f);
    
    %finds max likelyhood
    RSL2(i,1).pdfmax(:,1)=xi(index_to_max);
    RSL2(i,1).pdfmaxp(:,1)=f(index_to_max);    
    RSL2(i,1).pdfeu(:,1)=xi;
    RSL2(i,1).pdfeu(:,2)=f;
    %[phat,pci] = mle(xi,f);
    %RSL2(i,1).MLE_phat=phat;
    %RSL2(i,1).MLE_pci=pci;
           
    %B2(:,i)=dU2;       
    
    i=i+1;    
end
%B(end,:)=[];
%B2(end,:)=[];

% age-elevation linear regression
[p,s]=polyfit([RSL.age],[RSL.pdfmax],1); [f,delta] = polyval(p,[RSL.age],s);
R = corrcoef([RSL.age],[RSL.pdfmax]);
[pp,ss]=polyfit([RSL.age],[RSL.pr50],1); [ff,ddelta] = polyval(pp,[RSL.age],ss);

[p,s]=polyfit([RSL2.age],[RSL2.pdfmax],1); [f2,delta2] = polyval(p,[RSL2.age],s);
R2 = corrcoef([RSL2.age],[RSL2.pdfmax]);
[pp,ss]=polyfit([RSL2.age],[RSL2.pr50],1); [ff2,ddelta2] = polyval(pp,[RSL2.age],ss);

% outputs for LEM input and plotting figures
save plotdata/HRSL_data.mat

disp('..model ready..')
%% build RSL with 10 yr step for LEM model comparison
i=1; clear RSL B
age1=4;
age2=8.4; 
dage=0.01;

for age=age1:dage:age2
    mis1Zpred=[mis1.mis5eU] * age;
    dU=[mis1.Zrsl] - mis1Zpred;   
    RSL(i,1).age=age;
    RSL(i,1).eustatic=dU;
    pr = prctile(dU,[25 50 75]);
    RSL(i,1).pr25=pr(1);
    RSL(i,1).pr50=pr(2);
    RSL(i,1).pr75=pr(3);
    [f,xi] = ksdensity(dU);     
    [max_f, index_to_max] = max(f);    
    RSL(i,1).pdfmax(:,1)=xi(index_to_max);
    RSL(i,1).pdfmaxp(:,1)=f(index_to_max);    
    RSL(i,1).pdfeu(:,1)=xi;
    RSL(i,1).pdfeu(:,2)=f;           
    %B(:,i)=dU;
    i=i+1;    
end

save('plotdata/MIS1_RSL_10yr.mat','RSL')
disp('..HR model ready..')

%%
age1=6.2;
age2=7.4; 
RSL([RSL.age]<age1)=[];
RSL([RSL.age]>age2)=[];


%% average pdfs for GIA models peak ages
x_min = 0; x_max = 7;
n_points = 100;
x_eval = linspace(x_min, x_max, n_points);

i=1; 
age1=6.2;
age2=7.4; 
dage=0.1;

for age=age1:dage:age2
    mis1Zpred=[mis1.mis5eU] * age;
    dU=[mis1.Zrsl] - mis1Zpred;   
    [pdf_vals, x_vals] = ksdensity(dU, x_eval);
    pdfs(i,:) = pdf_vals;

    %dU2=[mis1.Zrsl_pre] - mis1Zpred;   
    %[pdf_vals2, x_vals2] = ksdensity(dU2, x_eval);    
    %pdfs2(i,:) = pdf_vals2;
    i=i+1;
end
    
%% Compute average PDF
HRSLp.x_eval=x_eval;
HRSLp.avg_pdf = mean(pdfs, 1);
HRSLp.avg_pdf2 = mean(pdfs2, 1);

% Compute confidence intervals (standard error)
age=age1:dage:age2;
HRSLp.n_datasets=numel(age);
HRSLp.std_pdf = std(pdfs, 0, 1);
HRSLp.std_pdf2 = std(pdfs2, 0, 1);
HRSLp.se_pdf = HRSLp.std_pdf / sqrt(numel(age));
HRSLp.se_pdf2 = HRSLp.std_pdf2 / sqrt(numel(age));
HRSLp.ci_lower = HRSLp.avg_pdf - 1.96 * HRSLp.se_pdf;
HRSLp.ci_upper = HRSLp.avg_pdf + 1.96 * HRSLp.se_pdf;
HRSLp.ci_lower2 = HRSLp.avg_pdf2 - 1.96 * HRSLp.se_pdf2;
HRSLp.ci_upper2 = HRSLp.avg_pdf2 + 1.96 * HRSLp.se_pdf2;

save('plotdata/HRSLp.mat','HRSLp')



