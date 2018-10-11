
include("./VGLib.jl");
using VGLib,HDF5,PyCall,PyPlot
@pyimport numpy as np

### Choose your cube here ###
dataname= "run-0-ppv.h5";
#---------------------------#

f       = h5open(dataname,"r");

ppv_d   = read( f, "PPV");    # PPV cube
v_ran   = read( f, "vrange"); # velocity range within the x-axis 
Q       = read( f, "Q");      # Q (Stokes Parameter)
U       = read( f, "U");      # U (Stokes Parameter)

I,C     = getIC_PPV(ppv_d,v_ran[1:end-1]); # Convert to the Intensity map and centroid map
cx,cy   = sobel_conv_2d(C);                # Find the x, y derivative of the centroid map
cna     = atan(cy./cx);                    # Calculating the magnetic using VCG
bs      = 72; 							   # Blocksize using in the averaging method 
Bcna,cns= sban2d(cos(cna),sin(cna),bs);    # Bloack averaging method 
pna     = .5*atan2(U,Q);                   # Calculating the magnetic direction from Pr
Bpna    = sban2d(cos(pna),sin(pna),bs)[1]; # Bloack averaging method 

imshow(I,cmap="binary")
Xd,Yd=np.meshgrid(div(bs,2):bs:size(I)[2],div(bs,2):bs:size(I)[1]);
quiver(Xd,Yd,cos(Bcna[:]),sin(Bcna[:]),headwidth=0,color="r",label="B field predicted by VGT");
quiver(Xd,Yd,-cos(Bcna[:]),-sin(Bcna[:]),headwidth=0,color="r")
quiver(Xd,Yd,cos(Bpna[:]),sin(Bpna[:]),headwidth=0,color="b",label="Projected B");
quiver(Xd,Yd,-cos(Bpna[:]),-sin(Bpna[:]),headwidth=0,color="b")
title("The comparsion between the VGT method and Projected B field direction")
axis("off");
legend(loc=3)
