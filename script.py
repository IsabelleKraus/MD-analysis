#!/usr/bin/env python
# coding: utf-8
# In[ ]:
#Compiled and modified by Dina Robaa

#RUN Script using: amber.python compiled-script.py

from __future__ import print_function
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) # avoid warning message to make the notebook nicer

#mkdir PLOTS to save generated files and Pics
import os

path = 'PLOTS'
# Check whether the specified path exists or not
isExist = os.path.exists(path)
if not isExist:
  # Create a new directory because it does not exist 
  os.makedirs(path)
  print("The new directory is created!")


#Import packages
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import colors
import pytraj as pt

#load trajectory
traj = pt.load('prod1.nc', top='z_solv.prmtop')				# If too many frames, the trajectory will fail to load -> memory error
#load input coordinates as reference structure for RMSD-calculation
ref = pt.load('z_solv.inpcrd', top=traj.top)
#align trajectory to input structure (ref)
traj.superpose(mask='@CA', ref=ref)

###RMSD for protein and ligand
##Here plotted vs Frame
#Use frame_indices=range(0, 20000, 50) for calculating only every 50th frame; see line 115
data_rmsd_pro = pt.rmsd(traj, ref=ref, mask="@CA")                  # Calculates RMSD for protein-backbone-CA; change mask for pro-residues! 
data_rmsd_lig = pt.rmsd(traj, ref=ref, mask=":LIG&!@H=", nofit=1)		# Calculates RMSD for ligand; change mask for residue-number 
frame_num=len(data_rmsd_pro)
m1_pro=max(data_rmsd_pro)
m1_lig=max(data_rmsd_lig)
#Save RMSD-data as csv
a=np.arange(1,frame_num+1,1)
data_frames=np.column_stack((a, data_rmsd_pro, data_rmsd_lig))
np.savetxt('PLOTS/RMSD_Frames.csv', data_frames, header='Frame\tprotein\tligand', fmt='%1.1f', delimiter='\t')
#RMSD-plot and save as RMSD_frames.png
plt.figure(0)                                                       # New Figure
plt.plot(data_rmsd_pro, label='Protein')	                          # Label legend
plt.plot(data_rmsd_lig, label='Ligand')		                          # Label legend
plt.xlabel('Frame', fontweight ='bold', fontsize=10)
plt.ylabel('RMSD (Angstrom)', fontweight ='bold', fontsize=10)		  # Label x and y axes 
plt.xlim(0,frame_num)				                                        # Range of x-axis
plt.xticks((np.arange(0,frame_num+1,50)), fontsize=10)              # Interval of ticks x-axis 
plt.yticks(fontsize=10)
plt.ylim(0,m1_lig+3)			                                          # range till max of lig-rmsd + 3 is used
plt.legend(fontsize=12)
plt.savefig('PLOTS/RMSD_frames.png', dpi=300)

###RMSD for protein and ligand
##Here plotted vs Time
T1=frame_num*0.02                                                       # Simulation Time. 0.02 has to be modified if ntwx changed. calculation 0.02: ntwx*dt/1000
Time=np.arange(0,T1,0.02)                                               # generate array for Time. Here 50 and 0.02 are chosen: (ntslim/ntwx/MD-Time)
data_pro_time=np.column_stack((Time,data_rmsd_pro))                     # data Time, rmsd-pro
data_lig_time=np.column_stack((Time,data_rmsd_lig))                     # data Time, rmsd-lig
data_time_all=np.column_stack((Time,data_rmsd_pro,data_rmsd_lig))       # data Time, rmsd-pro, rmsd-lig
#Save RMSD-data as csv
np.savetxt('PLOTS/RMSD_Time.csv', data_time_all, header='Time\tprotein\tligand', fmt='%1.1f', delimiter='\t')
#RMSD-plot vs Time and save as RMSD_Time.png
plt.figure(1)                                                     # New Figure
plt.plot(data_pro_time.T[0], data_pro_time.T[1], label='Protein')	# plot protein_RMSD and label as Protein
plt.plot(data_lig_time.T[0], data_lig_time.T[1], label='Ligand')	# plot protein_RMSD and label as Ligand
plt.xlabel('Time (ns)', fontweight ='bold', fontsize=10)
plt.ylabel('RMSD (Angstrom)', fontweight ='bold', fontsize=10)		# Label x and y axes
plt.xlim(0,T1)				                                            # range of x-axis
plt.xticks((np.arange(0,T1+1,1)), fontsize=10)                    # Interval of ticks x-axis
plt.yticks(fontsize=10)
plt.ylim(0,m1_lig+3)			                                        # range till max of lig-rmsd + 3 is used
plt.legend()
plt.savefig('PLOTS/RMSD_Time.png', dpi=300)


###RMSF-protein, save data and gerenerate RMSF-Plot:
rmsf_data = pt.rmsf(traj, mask="*&!(:WAT|:LIG,Na+,Cl-,K+,ZN)&!@H= byres")               
#rmsf_data.T[0]= rmsf_data.T[0] + 174       #Use this if protein sequence starts at 174 instead of at 1
residues=max(rmsf_data.T[0])
max_rmsf=max(rmsf_data.T[1])
np.savetxt('PLOTS/RMSF.csv', rmsf_data, fmt='%1.1f', delimiter='\t')
#Plot RMSF
plt.figure(2)                                                             
plt.plot(rmsf_data.T[0], rmsf_data.T[1], label='Protein')
plt.xlabel('Residue', fontweight ='bold', fontsize=10)
plt.ylabel('RMSF (Angstrom)', fontweight ='bold', fontsize=10)
plt.xlim(1,residues)
plt.ylim(0,max_rmsf+2)
plt.xticks((np.arange(0,residues+1,20)),rotation=50, fontsize=10)
plt.yticks(fontsize=10)
plt.legend()
plt.savefig('PLOTS/RMSF.png', dpi=300, bbox_inches='tight', pad_inches=0.2)

###RMSF-ligand-atoms,save data and gerenerate RMSF-Plot:
rmsf_data_lig = pt.rmsf(traj, mask=":LIG&!@H=") 
max_rmsf_lig=max(rmsf_data_lig.T[1])
np.savetxt('PLOTS/RMSF_lig.csv', rmsf_data_lig, fmt='%1.1f', delimiter='\t')
min_a_num=min(rmsf_data_lig.T[0])
max_a_num=max(rmsf_data_lig.T[0])

#Plot RMSF ligand
plt.figure(3)                                                       
plt.plot(rmsf_data_lig.T[0], rmsf_data_lig.T[1], label='Ligand')
plt.xlabel('Atom number', fontweight ='bold', fontsize=10)
plt.ylabel('RMSF (Angstrom)', fontweight ='bold', fontsize=10)
plt.xlim(min_a_num-1,max_a_num+1)
plt.xticks((np.arange(min_a_num,max_a_num+1)), rotation=50, fontsize=10)
plt.yticks(fontsize=10)
plt.ylim(0,max_rmsf_lig+2)
plt.legend()
plt.savefig('PLOTS/RMSF_lig.png', dpi=300, bbox_inches='tight', pad_inches=0.2)


###Compute pairwise RMSD of ligand (here for every 10th frame)
mat = pt.pairwise_rmsd(traj, mask=":LIG&!@H=", frame_indices=range(0, frame_num, 10))	#Every 100th frame is taken 
np.savetxt('PLOTS/Matrix.csv', mat, fmt='%1.1f', delimiter='\t')
plt.figure(4)
cmap = colors.ListedColormap(['royalblue', 'cornflowerblue', 'lightcoral', 'red'])
bounds = [0, 1, 2, 4, 6]
norm = colors.BoundaryNorm(bounds, cmap.N)
plt.imshow(mat, cmap=cmap, norm=norm)
plt.xlabel('Frame', fontweight ='bold', fontsize=10)
plt.ylabel('Frame', fontweight ='bold', fontsize=10)
plt.colorbar()
plt.savefig('PLOTS/Lig-RMSD-matrix.png', dpi=300)


### Hydrogen bonds Contacts
#Calculate Hbonds in trajectory
data= pt.search_hbonds(traj)
#Generate array with number of frames
string=pd.DataFrame({'frames'})
#Dataframe with found Hbonds
hbonds=pd.DataFrame(data.donor_acceptor, columns=['Hbond'])
#Generate Hbond matrix
array=pd.DataFrame(data.values[1:])
#Calculate number of frames
frames=array.shape[1]
#Concatenated dataframe with calculated Occupancy
df=pd.concat([hbonds, array], axis=1)
df.loc[:, 'Occupancy']=df.sum(axis=1)/frames*100
#Save all Hbonds found in trajectory to csv file, these include also intermolecular Hbonds
df.to_csv('PLOTS/Hbonds-all.csv', sep='\t')
#Generate dataframe with only lig-pro Hbonds
df1=df[df.iloc[:,0].str.contains('LIG')]
#Save Hbond matrix Lig-pro to csv file
df1.to_csv('PLOTS/Hbonds-Ligand.csv', sep='\t')
#Plot Hbond-matrix
column=list(df1.iloc[:,0])                                        # Create list with column 0
cm=cm = 1/2.54  # centimeters in inches
fig1, ax = plt.subplots(figsize=(25*cm,10*cm))   
                 # changed figure dimensions
ax.imshow(df1.iloc[:,1:], cmap='Blues',vmin=0, vmax=1, aspect='auto')
ax.set_yticks(np.arange(df1.shape[0]))
ax.set_yticklabels(list(df1.iloc[:,0]), fontsize=10)
plt.setp(ax.get_xticklabels(), ha="center", va="center")
plt.xlabel('Frame', fontsize=14, fontweight='bold')
ax.tick_params(top=False, bottom=True,labeltop=False, labelbottom=True)   # label of x-axis on bottom
fig1.tight_layout()
plt.savefig('PLOTS/H-bond-lig-pro.png', dpi=300)
#Save Hbond matrix Lig-pro to csv file
df1.to_csv('PLOTS/Hbonds-Ligand.csv', sep='\t')
#Save Hbond occupancy to csv
df1.to_csv('PLOTS/Hbonds-occupancy.csv', columns =['Hbond','Occupancy'], sep='\t')
Footer
