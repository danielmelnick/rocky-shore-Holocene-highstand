function LEM=GIALEMP(param,sl,UR,Z)

    param.uplift_rate=UR; %uplift rate m/ka
    
    %Create a LEM
    LEM.mod=TerraceM_LEM(sl,param);    
    LEM.UR=UR;    
    LEM.MIS1z=Z;         
    
    %find shoreline angle
    if numel([LEM.mod.x_mod])>0

        dxy=diffxy([LEM.mod.x_mod],[LEM.mod.z_mod]);
        ix=find(dxy==max(dxy)); ix=ix-1;   
    
        LEM.shx=[LEM.mod.x_mod(ix)];
        LEM.shz=[LEM.mod.z_mod(ix)];    
        LEM.dshz=LEM.shz-Z;    
    else
        disp('increase vx')
    end
end

