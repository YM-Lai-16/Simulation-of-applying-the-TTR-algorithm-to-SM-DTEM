# Simulation-of-applying-the-TTR-algorithm-to-SM-DTEM
Matlab code for implementing TwIST-based tomographic reconstruction (TTR) algorithm of compressed ultrafast tomographic imaging (CUTI) to the streak-mode dynamic transmission electron microscopy (SM-DTEM)

This script requires the package of the TwIST algorithm (TwIST_v2) available at 
http://www.lx.it.pt/~bioucas/TwIST/TwIST.htm

The customized functions "TTR_TV3D_conv.m", "TTR_TVphi.m", and "TTR_denoise_func.m" should be in the same folder of "TwIST.m"

The main file for running is "SMDTEM_TTR_Simulation.m".

Tested in Matlab R2020b
