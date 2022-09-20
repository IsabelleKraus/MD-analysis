#!/bin/bash
module load amber/22

cpptraj
parm z_solv.prmtop
trajin 01_Cosolvent/prod1.nc
trajin 02_Cosolvent/prod1.nc
trajin 03_Cosolvent/prod1.nc
trajin 04_Cosolvent/prod1.nc
trajin 05_Cosolvent/prod1.nc
reference 01_Cosolvent/z_solv.inpcrd
rms :1-55@CA reference
grid grid_EAM-N.dx 200 0.5 200 0.5 200 0.5 :EAM&@/N gridcenter 29.292695 33.157047 34.476436 pdb grid_EAM-N.pdb max 0.1   #use coordinates from bounds-inpcrd.dat center as gridcenter points
grid grid_EAM-O.dx 200 0.5 200 0.5 200 0.5 :EAM&@/O gridcenter 29.292695 33.157047 34.476436 pdb grid_EAM-O.pdb max 0.1   #0.1 means all occupancies over 10% 
run

