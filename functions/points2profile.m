function S = points2profile(line,points,fout)

% Project points along a line
%
% inputs
%
% line: shapefile in UTM
% point: shapefile in UTM
% fout: name of output file
%
% output 
% shapefile with field distance_along_profile in km
%
% D.Melnick 2008

S=shaperead(points); nim=numel(S);

Pf = shaperead(line); 
x1=Pf(1,1).X(1,1); x2=Pf(1,1).X(1,2); 
y1=Pf(1,1).Y(1,1); y2=Pf(1,1).Y(1,2);
A=[x1 y1]; B=[x2 y2];
AB = (B-A); 
AB_squared = dot(AB,AB); 

% get lat
[PLat,PLon] = utm18_2deg(Pf.X,Pf.Y);
AL=[PLon(1) PLat(1)]; 
BL=[PLon(2) PLat(2)];
ABL = (BL-AL); 
ABL_squared = dot(ABL,ABL); 

% project
stp=zeros(nim,4); 
for i=1:nim
    p=[S(i,1).X S(i,1).Y];
    Ap = (p-A); 
    t = dot(Ap,AB)/AB_squared; 
    p = A + t * AB;       
    stp(i,1)=p(1,1); 
    stp(i,2)=p(1,2);
    
    % calc dist
    stp(i,3)=hypot(p(1,1)-x1,p(1,2)-y1);    
    
    % calc projected Latitude
    [SLat,SLon] = utm18_2deg(S(i).X,S(i).Y);
    pL=[SLat SLon];
    Ap = (pL-AL); 
    t = dot(Ap,ABL)/ABL_squared; 
    p = AL + t * ABL;       
    %stp(i,4)=p(1,2) - Lat1; 
    
    % calc dist
    stp(i,4)=hypot(p(1,1)-PLon(1),p(1,2)-PLat(1));            
end

% add distance and projected latitude to output
for i=1:nim
    S(i).d_profile=stp(i,3);      
    S(i).d_profile_Lat=stp(i,4)-min(stp(:,4));
end

% save
shapewrite(S, fout)


