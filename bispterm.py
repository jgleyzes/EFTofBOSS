import numpy as np
#import matplotlib.pyplot as plt
import os.path as opa

#matplotlib inline
# THIS_PATH  =  opa.dirname(__file__)
thispath = '/exports/pierre/EFTofBOSS/input/GridsEFT/'
bispgrid = np.load(thispath+'TableBispconLightConeHectorPatchyWideHDvFFT_IR06.npy')

cosmologylist = np.load(thispath+'TablecoordconLightConeHectorwideh.npy')

params = cosmologylist[211276]
bispterms = bispgrid[211276]
kminbisp = 0.01
kmaxbisp = 0.07

Q1, Q2, Q3 = bispterms[:3]

masktriangle = (Q1>=kminbisp)&(Q1<=kmaxbisp)&(Q1<=Q2+Q3)&(Q1>=abs(Q2-Q3))&(Q2>=kminbisp)&(Q2<=kmaxbisp)&(Q3>=kminbisp)&(Q3<=kmaxbisp)

bval = np.array(['1.','b1','b2','b4','b1*b11','b1**2','b1*b2','b1*b4','b1**3','b1**2*b2','b1**2*b4','b8**2'])


for i in range(12):
    np.abs(bispterms[i+3][masktriangle]