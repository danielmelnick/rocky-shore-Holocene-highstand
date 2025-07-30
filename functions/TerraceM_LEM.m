function LEM=TerraceM_LEM(sl,param)

% Landscape Evolution Model for the generation of marine terraces after 
% Anderson et al. (2009) Basin Research, modified by Jara-Munoz et al
% (2019) Froniers in Earth Sciences for TerraceM Gui www.terracem.com

% INPUT:
% sl - Sea-level curve
% param - input structure with modeling setups
%
% OUTPUT:
% LEM - structure with model results

jj=1;

%parameters
tmax1=param.tmax;
tmin1=param.tmin;
dt=param.dt;
dx=param.dx;
slope_shelf=param.slope_shelf;
uplift_rate=param.uplift_rate;
edot0=param.initial_erosion;
kappa0=param.cliff_diffusion;
wave=param.wave_height;
vxx=param.vxx;
smth=param.smth;

sl=sl.ans;

% Can change the following variables: uplift rate (uprate) cliff erosion rate (edot0); and shelf slope (slope2)
if numel(uplift_rate)==1
    mode=1;
    uprate = uplift_rate/1000; % uplift rate in m/yr
else
    mode=2;
    uprate=max(uplift_rate)/1000;
    uprate666=flipud(uplift_rate./1000);
end

slope2 = tand(slope_shelf); 

tmax=tmax1*1000; %minimum  age 69000 / 83000
tmin=tmin1*1000; %maximum age 13000 /0

%%%%%%%%%%
M=sl(sl(:,1)<=tmax,:); 
M2=M(M(:,1)>=tmin,:); %adjust range for analysis in time window

%remove repeated
J=find(diff(M2(:,1))==0);
M2(J,:)=[];

xi = tmin:dt:tmax; %time steps
sea1=interp1(M2(:,1),M2(:,2),xi,'nearest')';
sea=smooth(-flipud(sea1)',smth);
time=xi;

%define Y axis range 
seamean=mean(sea);
seaamp=max(sea)-min(sea);
ttmax=max(time); 
ttmin=min(time);
jmax=((tmax-tmin)/dt);

% this is the depth to which one will feel waves significantly 
zstar=wave; %Wave height
wavedepth=3*zstar;
dissscale=dx*wavedepth;

slope = slope2; 
%zmin=(seamean-seaamp)-(tmax*uprate)-(vx/3) ;%-40
%zmax=(seamean+seaamp)+vx; %+200 %xmax=5000;
zmin=(seamean-seaamp)-(tmax*uprate)-seaamp-vxx ;%-40
zmax=(seamean+seaamp)+(tmax*uprate*1)+vxx;
xmax=(zmax-zmin)/slope;
x=1:dx:xmax; %x topography
imax=length(x);
z0=zmin+(slope*x);%y topography
z=z0;%y topography
xplot = x/1000; 
xmaxplot = max(xplot);

zsave=[x;z];
conto=1;

% ********** start the time clock ************* 
for j=1:jmax
% the uplift pattern, here taken to be uniform in x 
if mode==1
    uplift=uprate*dt*ones(1,length(x));    
else
    uplift=uprate666(j)*dt*ones(1,length(x));
end
surface=sea(j)*ones(1,length(x));
k=find(z>=surface);
z=z+uplift;

% find the intersection of the sea surface and the topography 
try
    xsea1(j)=x(k(1));
    tsea1(j)=(j)*dt;%time at time step
catch
    msgbox('Caution , low accuracy resut')
    TEXT=('Caution , low accuracy resut');
return
end
% dependence of the cliff erosion rate on the width of the shelf, now ignored
%Selected area for depth decay=depth of wave erosion
ggd=find(((sea(j)-wavedepth)<z)&(sea(j)>z)); %submarine part affected by waves

%boundaries of platform / experimental
length2_max=max(x(1,ggd));
lenght2_min=min(x(1,ggd));
length_total=length2_max-lenght2_min; %width of platform
edot1=edot0;
% the most important segments being closest to the instantaneous coastline
gg=find(((sea(j)-wavedepth)<z)&(sea(j)>z)); 
dissrate=exp((sea(j)-z(1,gg))/wavedepth)*dx;
dissipate=sum(dissrate);
edot=edot1*(dissscale/dissipate);
%edot=edot0;
ero(j)=edot;
ksea=k(1)+ceil(edot*dt/dx); %round reeplaced by ceil

try
  
shelf=z(k(1))+(slope2*(x-x(k(1))));
surface(1:(k(1)-1))=9999*ones(size(1:k(1)-1)); 
surface(ksea+1:imax)=9999*ones(size(ksea+1:imax)); 
z=min(z,surface);
xsea(j)=x(ksea+1); %-1

% now take into account erosion in this reach of
% nearshore area where water depths are < "wavedepth" 
ggnew=find(((sea(j)-wavedepth)<z)&(sea(j)>z)); 
nearshore=wavedepth/(x(max(ggnew))-x(ggnew(1))); 
znew=z(1,ggnew(1))+(nearshore*(x(1,ggnew)-x(1,ggnew(1)))); 
z(1,ggnew)=min(z(1,ggnew),znew);

%now diffuse the terrestrial portion of the space 
kappa=zeros(1,length(x));
kkk=find(x>xsea(j)); 
kappa1=kappa0*ones(1,length(kkk)); 
kappa(kkk(1):imax)=kappa1;

catch
msgbox('Error: increase the value Vx to extend the plot display') 
TEXT=('Info:  Error¡¡¡¡ Increase the value Vx to extend the plot display'); 
return
end 

dzdt = 4 * kappa .* ((del2(z)./(dx*dx))); 
z=z+(dzdt*dt);
tnow = (time(j)-tmax)/1000; 

xp1 = xplot;
xpatch = [0 xp1 xp1(imax:-1:1)];
zp1 = z;
zp2 = min(z0)*ones(1,length(x)); 
zpatch = [min(z0) zp1 zp2(imax:-1:1)];

LEM.prof(jj).t=tnow;
LEM.prof(jj).x=xpatch;
LEM.prof(jj).z=zpatch;
jj=jj+1;

conto=conto+1;
end

x_mod=x;
z_mod=z;
xsea1(end)=xsea1(end)+dx;%rectify xsea1 shoreline angle
idxs = arrayfun(@(x)find(x_mod==x,1),xsea1);
ysea1=z_mod(idxs); %maybe wrong

%LEM.xsea=xsea; 
LEM.xsea1=xsea1;
LEM.ysea1=ysea1;
LEM.x_mod=x;
LEM.z_mod=z;
LEM.sea_rec=max(tsea1)-tsea1;%tsea1
LEM.ero_mod=ero;
