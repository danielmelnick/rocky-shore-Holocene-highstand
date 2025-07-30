clear, clc
addpath(genpath(pwd))

% Load MIS1 shoreline database
load MIS1.mat

% Load SL
load D5i.mat 
load D6i.mat

% time 
t=D5i(:,1);
%% Compute a LEM for each Holocene shoreline using MIS5e UR for 5g and 6g GIA models
% uses parallel processing toolbox 

% set LEM input parameters
param.tmax=10;          % max time (ka)
param.tmin=0;           % min time (ka)
param.dt=10;            % bin size of time intervals (yrs)
param.dx=10;            % bin size of x axis (m)
param.slope_shelf=4;    % initial slope (deg)
param.uplift_rate=0.1;  % uplift rate m/ka
param.initial_erosion=0.01; % initial erosion rate E0 (m/yr)
param.cliff_diffusion=0.00; % cliff slope diffusion (m^2/yr) - set to zero for direct comparison with shoreline angle
param.wave_height=5;    % wave height (m)
param.vxx=500;          % extended x axis bound 
param.smth=0.0000001;   % smooth factor of sea-level curve (yr)

% ICE-6g
tic
for k=1:12    
    for i=1:numel(t)
        sl1(i,1)=t(i)*1e3;
        sl1(i,2)=-D6i(i,k+1);
        sl1(i,3)=1;     
    end
    sl.ans=sl1;
    % Run LEM for Holocene sites with MIS5e uplift rate
    parfor i=1:numel(mis1) 
        LEM6(i)=GIALEMP(param,sl,mis1(i).mis5eU,mis1(i).Zrsl); 
    end        
    GIAM6(k).mod=k;
    GIAM6(k).dshz=mean([LEM6.dshz]);
    GIAM6(k).sshz=std([LEM6.dshz]);
    GIAM6(k).adshz=mean(abs([LEM6.dshz]));
    GIAM6(k).sl=sl;
    %GIAM6(k).LEM6=LEM6;
    clear sl sl1    
    fprintf('...run %u/12 finished...\n',k)
end
toc

% ICE-5g
tic
for k=1:12
    for i=1:numel(t)
        sl2(i,1)=t(i)*1e3;
        sl2(i,2)=-D5i(i,k+1);
        sl2(i,3)=1;
    end
    sl.ans=sl2;
    % Run LEM for Holocene sites with MIS5e uplift rate
    parfor i=1:numel(mis1)        
        LEM5(i)=GIALEMP(param,sl,mis1(i).mis5eU,mis1(i).Zrsl);
    end 
    GIAM5(k).mod=k;
    GIAM5(k).dshz=mean([LEM5.dshz]);
    GIAM5(k).sshz=std([LEM5.dshz]);
    GIAM5(k).adshz=mean(abs([LEM5.dshz]));
    GIAM5(k).sl=sl;
    %GIAM5(k).LEM5=LEM5;

    clear sl sl2
    fprintf('...run %u/12 finished...\n',k)
end
toc

% add GIA model parameters to LEM outputs
for j=1:4
    GIAM6(j).V=j;
    GIAM6(j+4).V=j;
    GIAM6(j+8).V=j;
    GIAM6(j).L=71;
    GIAM6(j+4).L=96;
    GIAM6(j+8).L=120;
    GIAM5(j).L=71;
    GIAM5(j+4).L=96;
    GIAM5(j+8).L=120;        
    GIAM5(j).V=j;
    GIAM5(j+4).V=j;
    GIAM5(j+8).V=j;
end

% save outputs
save('indata/GIAM5_HRSL_LEM_postM.mat','GIAM5','-v7.3')
save('indata/GIAM6_HRSL_LEM_postM.mat','GIAM6','-v7.3')

%%
load GIAM5_HRSL_LEM_postM.mat
load GIAM6_HRSL_LEM_postM.mat
%
%% get best-fitting LEM model
a=[GIAM5.dshz];
b5=find(a==min(a));
a=[GIAM6.dshz];
b6=find(a==min(a));

% set LEM input parameters, see TerraceM for details
param.tmax=10;          % max time (ka)
param.tmin=0;           % min time (ka)
param.dt=10;            % bin size of time intervals (yrs)
param.dx=10;            % bin size of x axis (m)
param.slope_shelf=4;    % initial slope (deg)
param.uplift_rate=0.1;  % uplift rate m/ka
param.initial_erosion=0.01; % initial erosion rate E0 (m/yr)
param.cliff_diffusion=0.00; % cliff slope diffusion (m^2/yr) - set to zero for direct comparison with shoreline angle
param.wave_height=5;    % wave height (m)
param.vxx=500;          % extended x axis bound 
param.smth=0.0000001;   % smooth factor of sea-level curve (yr)

% Use SL of best-fit model
tic
for i=1:numel(t)
    sl5(i,1)=t(i)*1e3;
    sl5(i,2)=-D5i(i,b5+1);
    sl5(i,3)=1;
end
sl.ans=sl5;
% create set of LEM models with different URs using best GIA curve
parfor i=1:numel(mis1)            
    LEMB(i)=GIALEMP(param,sl,mis1(i).mis5eU,mis1(i).Zrsl);
end 
k=1;
mod(k).sl='5g';
mod(k).LEM=LEMB;
mod(k).mediands=median([LEMB.dshz]);
mod(k).meands=mean([LEMB.dshz]);
mod(k).std=std([LEMB.dshz]);

% 6g
for i=1:numel(t)
    sl6(i,1)=t(i)*1e3;    
    sl6(i,2)=-D6i(i,b6+1);
    sl6(i,3)=1;
end
sl.ans=sl6;
parfor i=1:numel(mis1)        
    LEMB(i)=GIALEMP(param,sl,mis1(i).mis5eU,mis1(i).Zrsl);    
end 
k=2;
mod(k).sl='6g';
mod(k).LEM=LEMB;
mod(k).mediands=median([LEMB.dshz]);
mod(k).meands=mean([LEMB.dshz]);
mod(k).std=std([LEMB.dshz]);
toc

%% save results
mod(1).LEM = rmfield(mod(1).LEM, 'mod');
mod(2).LEM = rmfield(mod(2).LEM, 'mod');
save('plotdata/HRSL_LEM_Best_postM.mat','mod','-v7.3')

