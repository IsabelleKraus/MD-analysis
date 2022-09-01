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
import re

#function to assign chain or ligand: (change according to residue numbers and names)
def assign_Chain(res_number):
    if res_number <= 55:
        result = 'A'
    elif res_number >= 58 and res_number <= 107:
        result = 'B'
    elif res_number >= 110 and res_number <= 162:
        result = 'C'
    elif res_number >= 165 and res_number <= 219:
        result = 'D'
    elif res_number >= 222 and res_number <= 226:
        result = 'LIG1'
    elif res_number >= 227 and res_number <= 230:
        result = 'LIG2'
    return result


#load trajectory
traj = pt.load('prod1.nc', top='z_solv.prmtop')				# If too many frames, the trajectory will fail to load -> memory error
#load input coordinates as reference structure for RMSD-calculation
ref = pt.load('z_solv.inpcrd', top=traj.top)
#Modify Trajectory 
#align trajectory to input protein structure (ref)
traj.superpose(mask=':1-221@CA', ref=ref)                           # changed with :1-221, to ignore ligand (peptide) atoms!

###############add: traj2 for defined superposing ####################################

###RMSD for protein and ligand
##Here plotted vs Frame
#Use frame_indices=range(0, 20000, 50) for calculating only every 50th frame; see line 115
data_rmsd_pro = pt.rmsd(traj, ref=ref, mask=":1-221@CA")            # Calculates RMSD for protein-backbone-CA; change mask for pro-residues! 
                                                                    # changed, see line 37
  
                                                                            # changed for 2 ligands and residue-numbers: 
data_rmsd_lig1 = pt.rmsd(traj, ref=ref, mask=":222-226&!@H=", nofit=True)		# Calculates RMSD for ligand; change mask for residue-number 
data_rmsd_lig2 = pt.rmsd(traj, ref=ref, mask=":227-230&!@H=", nofit=True)

frame_num=len(data_rmsd_pro)
m1_pro=max(data_rmsd_pro)
m1_lig1=max(data_rmsd_lig1)                                         # changed for 2 ligands
m1_lig2=max(data_rmsd_lig2)
m1_lig1_lig2=max(m1_lig1, m1_lig2)

plt.rcParams['figure.dpi'] = 600                                    # dots per inches: determines resolution when graph is plotted,
                                                                    # can be changed according to figure
#Save RMSD-data as csv
a=np.arange(1,frame_num+1,1)
data_frames=np.column_stack((a, data_rmsd_pro, data_rmsd_lig1, data_rmsd_lig2))     #changed
np.savetxt('PLOTS/RMSD_Frames.csv', data_frames, header='Frame\tprotein\tligand1\tligand2', fmt='%1.1f', delimiter='\t')
#RMSD-plot and save as RMSD_frames.png
plt.figure(0)                                                       # New Figure
plt.plot(data_rmsd_pro, label='Protein')	                          # Label legend
plt.plot(data_rmsd_lig1, label='Ligand 1')		                      # Label legend
plt.plot(data_rmsd_lig2, label='Ligand 2')                          # added for 2nd ligand, may change into sepparate plot
plt.xlabel('Frame', fontweight ='bold', fontsize=10)
plt.ylabel('RMSD (Angstrom)', fontweight ='bold', fontsize=10)		  # Label x and y axes 
plt.xlim(0,frame_num)				                                        # Range of x-axis
plt.xticks((np.arange(0,frame_num+1,500)), fontsize=10)              # Interval of ticks x-axis 
plt.yticks(fontsize=10)
plt.ylim(0,m1_lig1_lig2+3)			                                    # range till max of lig-rmsd + 3 is used #changed
plt.legend(fontsize=10)
plt.savefig('PLOTS/RMSD_frames.png', dpi=300)

###RMSD for protein and ligand
##Here plotted vs Time
T1=frame_num*0.1                                                       # Simulation Time. 0.02 has to be modified if ntwx changed. calculation 0.02: ntwx*dt/1000
Time=np.arange(0,T1,0.1)                                               # generate array for Time. Here 50 and 0.02 are chosen: (ntslim/ntwx/MD-Time)
data_pro_time=np.column_stack((Time,data_rmsd_pro))                     # data Time, rmsd-pro
data_lig1_time=np.column_stack((Time,data_rmsd_lig1))                   # data Time, rmsd-lig1                   # changed for 2 ligands
data_lig2_time=np.column_stack((Time,data_rmsd_lig2))
data_time_all=np.column_stack((Time,data_rmsd_pro,data_rmsd_lig1,data_rmsd_lig2))       # data Time, rmsd-pro, rmsd-lig
#Save RMSD-data as csv
np.savetxt('PLOTS/RMSD_Time.csv', data_time_all, header='Time\tprotein\tligand', fmt='%1.1f', delimiter='\t')
#RMSD-plot vs Time and save as RMSD_Time.png
plt.figure(1)                                                         # New Figure
plt.plot(data_pro_time.T[0], data_pro_time.T[1], label='Protein')	    # plot protein_RMSD and label as Protein
plt.plot(data_lig1_time.T[0], data_lig1_time.T[1], label='Ligand 1')	# plot protein_RMSD and label as Ligand     # changed for 2 ligands
plt.plot(data_lig2_time.T[0], data_lig2_time.T[1], label='Ligand 2')
plt.xlabel('Time (ns)', fontweight ='bold', fontsize=10)
plt.ylabel('RMSD (Angstrom)', fontweight ='bold', fontsize=10)		# Label x and y axes
plt.xlim(0,T1)				                                            # range of x-axis
plt.xticks((np.arange(0,T1+1,50)), fontsize=10)                    # Interval of ticks x-axis
plt.yticks(fontsize=10)
plt.ylim(0,m1_lig1_lig2+3)			                                        # range till max of lig-rmsd + 3 is used
plt.legend()
plt.savefig('PLOTS/RMSD_Time.png', dpi=300)

###RMSF-protein, save data and gerenerate RMSF-Plot:
rmsf_data = pt.rmsf(traj, mask=":1-221@CA byres")               
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
plt.axvspan(0, 57, alpha=0.1, color='green')                                # to color parts of plot according to protein chain 
plt.axvspan(58, 107, alpha=0.1, color='orange')
plt.axvspan(110, 162, alpha=0.1, color='green')
plt.axvspan(165, 219, alpha=0.1, color='orange')
plt.savefig('PLOTS/RMSF.png', dpi=300, bbox_inches='tight', pad_inches=0.2)




###RMSF-ligand1-atoms,save data and gerenerate RMSF-Plot:
rmsf_data_lig = pt.rmsf(traj, mask=":222-226&!@H= byres") 
max_rmsf_lig=max(rmsf_data_lig.T[1])
np.savetxt('PLOTS/RMSF_lig_1.csv', rmsf_data_lig, fmt='%1.1f', delimiter='\t')
min_a_num=min(rmsf_data_lig.T[0])
max_a_num=max(rmsf_data_lig.T[0])

#Plot RMSF ligand1
plt.figure(3)                                                       
plt.plot(rmsf_data_lig.T[0], rmsf_data_lig.T[1], label='Ligand 1')
plt.xlabel('Residue number', fontweight ='bold', fontsize=10)
plt.ylabel('RMSF (Angstrom)', fontweight ='bold', fontsize=10)
plt.xlim(min_a_num-1,max_a_num+1)
plt.xticks((np.arange(min_a_num,max_a_num+1)), rotation=50, fontsize=10)
plt.yticks(fontsize=10)
plt.ylim(0,max_rmsf_lig+2)
plt.legend()
plt.savefig('PLOTS/RMSF_lig_1.png', dpi=300, bbox_inches='tight', pad_inches=0.2)

###RMSF-ligand2-atoms,save data and gerenerate RMSF-Plot:
rmsf_data_lig_2 = pt.rmsf(traj, mask=":58-61&!@H= byres") 
max_rmsf_lig_2=max(rmsf_data_lig_2.T[1])
np.savetxt('PLOTS/RMSF_lig_2.csv', rmsf_data_lig_2, fmt='%1.1f', delimiter='\t')
min_a_num=min(rmsf_data_lig_2.T[0])
max_a_num=max(rmsf_data_lig_2.T[0])

#Plot RMSF ligand2
plt.figure(3)                                                       
plt.plot(rmsf_data_lig_2.T[0], rmsf_data_lig_2.T[1], label='Ligand 2')
plt.xlabel('Residue number', fontweight ='bold', fontsize=10)
plt.ylabel('RMSF (Angstrom)', fontweight ='bold', fontsize=10)
plt.xlim(min_a_num-1,max_a_num+1)
plt.xticks((np.arange(min_a_num,max_a_num+1)), rotation=50, fontsize=10)
plt.yticks(fontsize=10)
plt.ylim(0,max_rmsf_lig_2+2)
plt.legend()
plt.savefig('PLOTS/RMSF_lig_2.png', dpi=300, bbox_inches='tight', pad_inches=0.2)


### Hydrogen bonds Contacts
#Calculate Hbonds in trajectory
data= pt.search_hbonds(traj)
#Generate array with number of frames
string=pd.DataFrame({'frames'})
#Dataframe with found Hbonds
hbonds=pd.DataFrame(data.donor_acceptor, columns=['Hbond'])
#print(hbonds)
#Generate Hbond matrix
array=pd.DataFrame(data.values[1:])
#Calculate number of frames
frames=array.shape[1]
#print(frames)
#Concatenated dataframe with calculated Occupancy
df=pd.concat([hbonds, array], axis=1)
df.loc[:, 'Occupancy']=df.sum(axis=1)/frames*100
print(df)
#Save all Hbonds found in trajectory to csv file, these include also intermolecular Hbonds
df.to_csv('PLOTS/Hbonds-all.csv', sep='\t')

#Generate dataframe with only lig1-pro Hbonds
df1=df[df.iloc[:,0].str.contains('222|223|224|225|226')]
#Generate dataframe with lig-pro H-bond Occupancy over 10%
df1=df1.loc[df1['Occupancy']>=10]
df1=df1.sort_values(by=['Hbond'])
#Generate 2 columns with Hbond involved residues 1 and 2
df1=df1.assign(Res1=df1['Hbond'].str[:6])
df1['Res2']=df1['Hbond'].str.split('-').str[1].str[:6]
#Sort dataframe after first residue and then occupancy
df1=df1.sort_values(by=['Res1','Occupancy'])
#extract number of residue and add it to dataframe
df1['ChainRes1']=df1['Res1'].str.extract(r'([0-9]+)')
df1['ChainRes2']=df1['Res2'].str.extract(r'([0-9]+)')
#cast residue number to integer
df1['ChainRes1']=pd.to_numeric(df1['ChainRes1'])
df1['ChainRes2']=pd.to_numeric(df1['ChainRes2'])
#change residue number into assigned chain or ligand 
df1['ChainRes1'] = df1['ChainRes1'].apply(lambda x: assign_Chain(x))
df1['ChainRes2'] = df1['ChainRes2'].apply(lambda x: assign_Chain(x))
#add assignd chain/ligand to name of Hbond
df1['Hbond']=df1['Hbond']+' ('+ df1['ChainRes1'] + ',' + df1['ChainRes2'] + ')' 
print(df1)
df1.to_csv('PLOTS/Hbonds-Lig1_extended.csv', columns=['Hbond','Occupancy','Res1','Res2','ChainRes1','ChainRes2'], sep='\t')
#delete 2 residue columns used for sorting and chain assignments
df1=df1.drop(columns=['Res1','Res2','ChainRes1','ChainRes2'])
#print(df1)

#Save Hbond matrix Lig-pro to csv file
df1.to_csv('PLOTS/Hbonds-Ligand_1.csv', sep='\t')
#Plot Hbond-matrix
plt.rcParams['figure.dpi'] = 1000
column=list(df1.iloc[:,0])                                        # Create list with column 0
cm=cm = 1/2.54  # centimeters in inches
fig1, ax = plt.subplots(figsize=(50*cm,30*cm))   
                 # changed figure dimensions
ax.imshow(df1.iloc[:,1::5], cmap='Blues',vmin=0, vmax=1, aspect='auto')
ax.set_yticks(np.arange(df1.shape[0]))
ax.set_yticklabels(list(df1.iloc[:,0]), fontsize=10)
plt.setp(ax.get_xticklabels(), ha="center", va="center")
plt.xlabel('Frame', fontsize=14, fontweight='bold')
ax.tick_params(top=False, bottom=True,labeltop=False, labelbottom=True)   # label of x-axis on bottom

_label_list = [1000,2000,3000,4000,5000]
ax.set_xticks([0,200,400,600,800,1000])
ax.set_xticklabels(x_label_list)

fig1.tight_layout()
plt.savefig('PLOTS/H-bond-lig_1-pro.png', dpi=300)
#Save Hbond matrix Lig-pro to csv file
df1.to_csv('PLOTS/Hbonds-Ligand_1.csv', sep='\t')
#Save Hbond occupancy to csv
df1.to_csv('PLOTS/Hbonds-occupancy_Lig_1.csv', columns =['Hbond','Occupancy'], sep='\t')

#Generate dataframe with only lig_2-pro Hbonds
df2=df[df.iloc[:,0].str.contains('227|228|229|230')]
df2=df2.loc[df2['Occupancy']>=10]                                         # generate dataframe with lig-pro H-bond Occupancy over 10%
df2=df2.sort_values(by=['Hbond'])                                         # generate 2 columns with Hbond involved residues 1 and 2
df2=df2.assign(Res1=df2['Hbond'].str[:6])
df2['Res2']=df2['Hbond'].str.split('-').str[1].str[:6]
df2=df2.sort_values(by=['Res1','Occupancy'])                              # sort dataframe after first residue and occupancy
df2['ChainRes1']=df2['Res1'].str.extract(r'([0-9]+)')                     # extract number of residue and add it to dataframe
df2['ChainRes2']=df2['Res2'].str.extract(r'([0-9]+)')
df2['ChainRes1']=pd.to_numeric(df2['ChainRes1'])                          # cast residue number to integer
df2['ChainRes2']=pd.to_numeric(df2['ChainRes2'])
df2['ChainRes1'] = df2['ChainRes1'].apply(lambda x: assign_Chain(x))      # change residue numbers in columns to assigned chain
df2['ChainRes2'] = df2['ChainRes2'].apply(lambda x: assign_Chain(x))
df2['Hbond']=df2['Hbond']+' ('+ df2['ChainRes1'] + ',' + df2['ChainRes2'] + ')' 
print(df2)
df2.to_csv('PLOTS/Hbonds-Lig2_extended.csv', columns=['Hbond','Occupancy','Res1','Res2','ChainRes1','ChainRes2'], sep='\t')   #save extended dataframe
df2=df2.drop(columns=['Res1','Res2','ChainRes1','ChainRes2'])             # delete 4 residue columns used for sorting and chain assignments

#Save Hbond matrix Lig2-pro to csv file
df2.to_csv('PLOTS/Hbonds-Ligand_2.csv', sep='\t')
#Plot Hbond-matrix
column=list(df2.iloc[:,0])                                        # Create list with column 0
cm=cm = 1/2.54  # centimeters in inches
fig1, ax = plt.subplots(figsize=(50*cm,30*cm))   
                 # changed figure dimensions
ax.imshow(df2.iloc[:,1:], cmap='Blues',vmin=0, vmax=1, aspect='auto')
ax.set_yticks(np.arange(df2.shape[0]))
ax.set_yticklabels(list(df2.iloc[:,0]), fontsize=10)
plt.setp(ax.get_xticklabels(), ha="center", va="center")
plt.xlabel('Frame', fontsize=14, fontweight='bold')
ax.tick_params(top=False, bottom=True,labeltop=False, labelbottom=True)   # label of x-axis on bottom
fig1.tight_layout()
plt.savefig('PLOTS/H-bond-lig_2-pro.png', dpi=600)
#Save Hbond matrix Lig-pro to csv file
df2.to_csv('PLOTS/Hbonds-Ligand_2.csv', sep='\t')
#Save Hbond occupancy to csv
df2.to_csv('PLOTS/Hbonds-occupancy_Lig_2.csv', columns =['Hbond','Occupancy'], sep='\t')


###RMSD for protein chains
##Here plotted vs Frame
#Use frame_indices=range(0, 20000, 50) for calculating only every 50th frame; see line 115
data_rmsd_proA = pt.rmsd(traj, ref=ref, mask=":1-57@CA", nofit=False)# Calculates RMSD for protein-backbone-CA; change mask for pro-residues! 
data_rmsd_proB = pt.rmsd(traj, ref=ref, mask=":58-107@CA", nofit=False)
data_rmsd_proC = pt.rmsd(traj, ref=ref, mask=":110-162@CA", nofit=False)
data_rmsd_proD = pt.rmsd(traj, ref=ref, mask=":165-219@CA", nofit=False)
  
                                                                          # changed for 2 ligands and residue-numbers: 
#data_rmsd_lig1 = pt.rmsd(traj, ref=ref, mask=":222-226&!@H=", nofit=False)		# Calculates RMSD for ligand; change mask for residue-number 
#data_rmsd_lig2 = pt.rmsd(traj, ref=ref, mask=":227-230&!@H=", nofit=False)

#pl.rcParams['figure.dpi'] = 300 # reset size of plot window

frame_num=len(data_rmsd_proA)
m1_pro=max(data_rmsd_proC)
m1_lig1=max(data_rmsd_lig1)                                         # changed for 2 ligands
m1_lig2=max(data_rmsd_lig2)
m1_lig1_lig2=max(m1_lig1, m1_lig2)

#Save RMSD-data as csv
a=np.arange(1,frame_num+1,1)
data_frames=np.column_stack((a, data_rmsd_proA, data_rmsd_proB, data_rmsd_proC, data_rmsd_proD))     #changed
np.savetxt('PLOTS/RMSD_Frames_ch.csv', data_frames, header='Frame\tprotein_chainA\tprotein_chainB\tprotein_chainC\tprotein_chainD\tligand1\tligand2', fmt='%1.1f', delimiter='\t')
#RMSD-plot and save as RMSD_frames.png
plt.figure(0)                                                       # New Figure
plt.plot(data_rmsd_proA, label='Protein Chain A')
plt.plot(data_rmsd_proB, label='Protein Chain B')
plt.plot(data_rmsd_proC, label='Protein Chain C')
plt.plot(data_rmsd_proD, label='Protein Chain D')
plt.xlabel('Frame', fontweight ='bold', fontsize=10)
plt.ylabel('RMSD (Angstrom)', fontweight ='bold', fontsize=10)		  # Label x and y axes 
plt.xlim(0,frame_num)				                                        # Range of x-axis
plt.xticks((np.arange(0,frame_num+1,500)), fontsize=10)              # Interval of ticks x-axis 
plt.yticks(fontsize=10)
plt.ylim(0,m1_pro+3)			                                    # range till max of lig-rmsd + 3 is used #changed
plt.legend(fontsize=12)
plt.savefig('PLOTS/RMSD_frames_ch.png', dpi=300)

###RMSD for protein chains
##Here plotted vs Time
T1=frame_num*0.1                                                       # Simulation Time. 0.02 has to be modified if ntwx changed. calculation 0.02: ntwx*dt/1000
Time=np.arange(0,T1,0.1)                                               # generate array for Time. Here 50 and 0.02 are chosen: (ntslim/ntwx/MD-Time)
data_pro_timeA=np.column_stack((Time,data_rmsd_proA))                   # data Time, rmsd-pro
data_pro_timeB=np.column_stack((Time,data_rmsd_proB))
data_pro_timeC=np.column_stack((Time,data_rmsd_proC))
data_pro_timeD=np.column_stack((Time,data_rmsd_proD))

data_lig1_time=np.column_stack((Time,data_rmsd_lig1))                   # data Time, rmsd-lig1                   # changed for 2 ligands
data_lig2_time=np.column_stack((Time,data_rmsd_lig2))
data_time_all=np.column_stack((Time,data_rmsd_proA,data_rmsd_proB,data_rmsd_proC,data_rmsd_proD,data_rmsd_lig1,data_rmsd_lig2))       # data Time, rmsd-pro, rmsd-lig
#Save RMSD-data as csv
np.savetxt('PLOTS/RMSD_Time.csv', data_time_all, header='Time\tprotein_chainA\tprotein_chainB\tprotein_chainC\tprotein_chainD\tligand1\tligand2', fmt='%1.1f', delimiter='\t')
#RMSD-plot vs Time and save as RMSD_Time.png
plt.figure(1)                                                         # New Figure
plt.plot(data_pro_timeA.T[0], data_pro_timeA.T[1], label='Protein Chain A')	    # plot protein_RMSD and label as Protein
plt.plot(data_pro_timeB.T[0], data_pro_timeB.T[1], label='Protein Chain B')
plt.plot(data_pro_timeC.T[0], data_pro_timeC.T[1], label='Protein Chain C')
plt.plot(data_pro_timeD.T[0], data_pro_timeD.T[1], label='Protein Chain D')
#plt.plot(data_lig1_time.T[0], data_lig1_time.T[1], label='Ligand 1')	# plot protein_RMSD and label as Ligand     # changed for 2 ligands
#plt.plot(data_lig2_time.T[0], data_lig2_time.T[1], label='Ligand 2')
plt.xlabel('Time (ns)', fontweight ='bold', fontsize=10)
plt.ylabel('RMSD (Angstrom)', fontweight ='bold', fontsize=10)		# Label x and y axes
plt.xlim(0,T1)				                                            # range of x-axis
plt.xticks((np.arange(0,T1+1,50)), fontsize=10)                    # Interval of ticks x-axis
plt.yticks(fontsize=10)
plt.ylim(0,m1_pro+3)			                                        # range till max of lig-rmsd + 3 is used
plt.legend()
plt.savefig('PLOTS/RMSD_Time_ch.png', dpi=300)


#pl.rcParams['figure.dpi'] = 300 # reset size of plot window


###RMSD for protein chain B,D and lig1
##Here plotted vs Frame
#Use frame_indices=range(0, 20000, 50) for calculating only every 50th frame; see line 115
#nofit=False for protein to see degree of intern movements
data_rmsd_pro = pt.rmsd(traj, ref=ref, mask=":58-107,165-219@CA", nofit=False)# Calculates RMSD for protein-backbone-CA; change mask for pro-residues! 
#data_rmsd_proB = pt.rmsd(traj, ref=ref, mask=":58-107@CA", nofit=True)
#data_rmsd_proC = pt.rmsd(traj, ref=ref, mask=":110-162@CA", nofit=True)
#data_rmsd_proD = pt.rmsd(traj, ref=ref, mask=":165-219@CA", nofit=True)
  
                                                                          # changed for 2 ligands and residue-numbers: 
data_rmsd_lig1 = pt.rmsd(traj, ref=ref, mask=":222-226&!@H=", nofit=True)		# Calculates RMSD for ligand; change mask for residue-number 
data_rmsd_lig2 = pt.rmsd(traj, ref=ref, mask=":227-230&!@H=", nofit=True)

frame_num=len(data_rmsd_pro)
m1_pro=max(data_rmsd_pro)
m1_lig1=max(data_rmsd_lig1)                                         # changed for 2 ligands
m1_lig2=max(data_rmsd_lig2)
m1_lig1_lig2=max(m1_lig1, m1_lig2)



#Save RMSD-data as csv
a=np.arange(1,frame_num+1,1)
data_frames=np.column_stack((a, data_rmsd_pro , data_rmsd_lig1))     #changed
np.savetxt('PLOTS/RMSD_Frames_chBD.csv', data_frames, header='Frame\tprotein_chainA\tprotein_chainB\tprotein_chainC\tprotein_chainD\tligand1\tligand2', fmt='%1.1f', delimiter='\t')
#RMSD-plot and save as RMSD_frames.png
plt.figure(0)                                                       # New Figure
plt.plot(data_rmsd_pro, label='Protein Chain B and D')
#plt.plot(data_rmsd_proB, label='Protein Chain B')
#plt.plot(data_rmsd_proC, label='Protein Chain C')
#plt.plot(data_rmsd_proD, label='Protein Chain D')

plt.plot(data_rmsd_lig1, label='Ligand 1')		                    # Label legend
#plt.plot(data_rmsd_lig2, label='Ligand 2')                          # added for 2nd ligand, may change into sepparate plot
plt.xlabel('Frame', fontweight ='bold', fontsize=10)
plt.ylabel('RMSD (Angstrom)', fontweight ='bold', fontsize=10)		  # Label x and y axes 
plt.xlim(0,frame_num)				                                        # Range of x-axis
plt.xticks((np.arange(0,frame_num+1,500)), fontsize=10)              # Interval of ticks x-axis 
plt.yticks(fontsize=10)
plt.ylim(0,m1_pro+3)			                                    # range till max of lig-rmsd + 3 is used #changed
plt.legend(fontsize=12)
plt.savefig('PLOTS/RMSD_frames_chBD.png', dpi=300)

###RMSD for protein chaind B, D and ligand1
##Here plotted vs Time
T1=frame_num*0.1                                                       # Simulation Time. 0.02 has to be modified if ntwx changed. calculation 0.02: ntwx*dt/1000
Time=np.arange(0,T1,0.1)                                               # generate array for Time. Here 50 and 0.02 are chosen: (ntslim/ntwx/MD-Time)
data_pro_timeA=np.column_stack((Time,data_rmsd_pro))                   # data Time, rmsd-pro
#data_pro_timeB=np.column_stack((Time,data_rmsd_proB))
#data_pro_timeC=np.column_stack((Time,data_rmsd_proC))
#data_pro_timeD=np.column_stack((Time,data_rmsd_proD))

data_lig1_time=np.column_stack((Time,data_rmsd_lig1))                   # data Time, rmsd-lig1                   # changed for 2 ligands
#data_lig2_time=np.column_stack((Time,data_rmsd_lig2))
data_time_all=np.column_stack((Time,data_rmsd_pro,data_rmsd_lig1))       # data Time, rmsd-pro, rmsd-lig
#Save RMSD-data as csv
np.savetxt('PLOTS/RMSD_Time_chBD.csv', data_time_all, header='Time\tprotein_chainBD\tligand1', fmt='%1.1f', delimiter='\t')
#RMSD-plot vs Time and save as RMSD_Time.png
plt.figure(1)                                                         # New Figure
plt.plot(data_pro_timeA.T[0], data_pro_timeA.T[1], label='Protein Chain B and D')	    # plot protein_RMSD and label as Protein
#plt.plot(data_pro_timeB.T[0], data_pro_timeB.T[1], label='Protein Chain B')
#plt.plot(data_pro_timeC.T[0], data_pro_timeC.T[1], label='Protein Chain C')
#plt.plot(data_pro_timeD.T[0], data_pro_timeD.T[1], label='Protein Chain D')

plt.plot(data_lig1_time.T[0], data_lig1_time.T[1], label='Ligand 1')	# plot protein_RMSD and label as Ligand     # changed for 2 ligands
#plt.plot(data_lig2_time.T[0], data_lig2_time.T[1], label='Ligand 2')
plt.xlabel('Time (ns)', fontweight ='bold', fontsize=10)
plt.ylabel('RMSD (Angstrom)', fontweight ='bold', fontsize=10)		# Label x and y axes
plt.xlim(0,T1)				                                            # range of x-axis
plt.xticks((np.arange(0,T1+1,50)), fontsize=10)                    # Interval of ticks x-axis
plt.yticks(fontsize=10)
plt.ylim(0,m1_pro+3)			                                        # range till max of lig-rmsd + 3 is used
plt.legend()
plt.savefig('PLOTS/RMSD_Time_chBD.png', dpi=300)


#For RMSF plot where spcific chains are superposed (aligned) --> tranlation and rotation not considered
traj2 = traj[:] #copy of original trajectory
traj2.superpose(mask=':85-107,165-219,222-226@CA') #align copied trajectory to initial structure of chain B, D and ligand 1

###RMSF-ligand-atoms with superposition,save data and generate RMSF-Plot:
rmsf_data_lig = pt.rmsf(traj2, mask=":222-226&!@H= byres") 
max_rmsf_lig=max(rmsf_data_lig.T[1])
np.savetxt('PLOTS/RMSF_lig_1_traj2.csv', rmsf_data_lig, fmt='%1.1f', delimiter='\t')
min_a_num=min(rmsf_data_lig.T[0])
max_a_num=max(rmsf_data_lig.T[0])

#Plot RMSF ligand
plt.figure(3)                                                       
plt.plot(rmsf_data_lig.T[0], rmsf_data_lig.T[1], label='Ligand 1')
plt.xlabel('Residue number', fontweight ='bold', fontsize=10)
plt.ylabel('RMSF (Angstrom)', fontweight ='bold', fontsize=10)
plt.xlim(min_a_num-1,max_a_num+1)
plt.xticks((np.arange(min_a_num,max_a_num+1)), rotation=50, fontsize=10)
plt.yticks(fontsize=10)
plt.ylim(0,max_rmsf_lig+2)
plt.legend()
plt.savefig('PLOTS/RMSF_lig_1_traj2.png', dpi=300, bbox_inches='tight', pad_inches=0.2)

#pl.rcParams['figure.dpi'] = 300 # reset size of plot window


###RMSD for protein chain A,C and lig2
##Here plotted vs Frame
#Use frame_indices=range(0, 20000, 50) for calculating only every 50th frame; see line 115
#nofit=False for protein to see degree of intern movements
data_rmsd_proAC = pt.rmsd(traj, ref=ref, mask=":110-162,165-219@CA", nofit=False)# Calculates RMSD for protein-backbone-CA; change mask for pro-residues! 
#data_rmsd_proB = pt.rmsd(traj, ref=ref, mask=":58-107@CA", nofit=True)
#data_rmsd_proC = pt.rmsd(traj, ref=ref, mask=":110-162@CA", nofit=True)
#data_rmsd_proD = pt.rmsd(traj, ref=ref, mask=":165-219@CA", nofit=True)
  
                                                                          # changed for 2 ligands and residue-numbers: 
data_rmsd_lig1 = pt.rmsd(traj, ref=ref, mask=":222-226&!@H=", nofit=True)		# Calculates RMSD for ligand; change mask for residue-number 
data_rmsd_lig2 = pt.rmsd(traj, ref=ref, mask=":227-230&!@H=", nofit=True)

frame_num=len(data_rmsd_proAC)
m1_pro=max(data_rmsd_proAC)
m1_lig1=max(data_rmsd_lig1)                                         # changed for 2 ligands
m1_lig2=max(data_rmsd_lig2)
m1_lig1_lig2=max(m1_lig1, m1_lig2)

#Save RMSD-data as csv
a=np.arange(1,frame_num+1,1)
data_frames=np.column_stack((a, data_rmsd_proAC , data_rmsd_lig2))     #changed
np.savetxt('PLOTS/RMSD_Frames_chAC.csv', data_frames, header='Frame\tprotein_chainA\tprotein_chainC\tligand2', fmt='%1.1f', delimiter='\t')
#RMSD-plot and save as RMSD_frames.png
plt.figure(0)                                                       # New Figure
plt.plot(data_rmsd_proAC, label='Protein Chain A and C')
#plt.plot(data_rmsd_proB, label='Protein Chain B')
#plt.plot(data_rmsd_proC, label='Protein Chain C')
#plt.plot(data_rmsd_proD, label='Protein Chain D')

plt.plot(data_rmsd_lig2, label='Ligand 2')		                    # Label legend
#plt.plot(data_rmsd_lig2, label='Ligand 2')                          # added for 2nd ligand, may change into sepparate plot
plt.xlabel('Frame', fontweight ='bold', fontsize=10)
plt.ylabel('RMSD (Angstrom)', fontweight ='bold', fontsize=10)		  # Label x and y axes 
plt.xlim(0,frame_num)				                                        # Range of x-axis
plt.xticks((np.arange(0,frame_num+1,500)), fontsize=10)              # Interval of ticks x-axis 
plt.yticks(fontsize=10)
plt.ylim(0,m1_lig2+3)			                                    # range till max of lig-rmsd + 3 is used #changed
plt.legend(fontsize=12)
plt.savefig('PLOTS/RMSD_frames_chAC.png', dpi=300)


###RMSD for protein chaind A, C and ligand2
##Here plotted vs Time
T1=frame_num*0.1                                                       # Simulation Time. 0.02 has to be modified if ntwx changed. calculation 0.02: ntwx*dt/1000
Time=np.arange(0,T1,0.1)                                               # generate array for Time. Here 50 and 0.02 are chosen: (ntslim/ntwx/MD-Time)
data_pro_timeA=np.column_stack((Time,data_rmsd_proAC))                   # data Time, rmsd-pro
#data_pro_timeB=np.column_stack((Time,data_rmsd_proB))
#data_pro_timeC=np.column_stack((Time,data_rmsd_proC))
#data_pro_timeD=np.column_stack((Time,data_rmsd_proD))

#data_lig1_time=np.column_stack((Time,data_rmsd_lig1))                   # data Time, rmsd-lig1                   # changed for 2 ligands
data_lig2_time=np.column_stack((Time,data_rmsd_lig2))
data_time_all=np.column_stack((Time,data_rmsd_proAC,data_rmsd_lig2))       # data Time, rmsd-pro, rmsd-lig
#Save RMSD-data as csv
np.savetxt('PLOTS/RMSD_Time_chAC.csv', data_time_all, header='Time\tprotein_chainAC\tligand1', fmt='%1.1f', delimiter='\t')
#RMSD-plot vs Time and save as RMSD_Time.png
plt.figure(1)                                                         # New Figure
plt.plot(data_pro_timeA.T[0], data_pro_timeA.T[1], label='Protein Chain A and C')	    # plot protein_RMSD and label as Protein
#plt.plot(data_pro_timeB.T[0], data_pro_timeB.T[1], label='Protein Chain B')
#plt.plot(data_pro_timeC.T[0], data_pro_timeC.T[1], label='Protein Chain C')
#plt.plot(data_pro_timeD.T[0], data_pro_timeD.T[1], label='Protein Chain D')

#plt.plot(data_lig1_time.T[0], data_lig1_time.T[1], label='Ligand 1')	# plot protein_RMSD and label as Ligand     # changed for 2 ligands
plt.plot(data_lig2_time.T[0], data_lig2_time.T[1], label='Ligand 2')
plt.xlabel('Time (ns)', fontweight ='bold', fontsize=10)
plt.ylabel('RMSD (Angstrom)', fontweight ='bold', fontsize=10)		# Label x and y axes
plt.xlim(0,T1)				                                            # range of x-axis
plt.xticks((np.arange(0,T1+1,50)), fontsize=10)                    # Interval of ticks x-axis
plt.yticks(fontsize=10)
plt.ylim(0,m1_lig2+3)			                                        # range till max of lig-rmsd + 3 is used
plt.legend()
plt.savefig('PLOTS/RMSD_Time_chAC.png', dpi=300)



#RMSF-Protein-Plot with RMSF for all chains and right residue numbers 

seq= "NMVHPNVICDGCNGPVVGTRYKCSVCPDYDLCSVCEGKGLHRGHTKLAFPSPF"    # aligned sequence (every chain just differs in length)
aa= list(seq)
res =[]

res_num=120                                                     # to concat residue number and residue name
for j in aa:
    j = j + " " + str(res_num)
    res.append(j)
    res_num += 1

traj.superpose(mask=':1-221@CA', ref=ref)   
rmsf_dataA = pt.rmsf(traj, mask=":1-55@CA byres")
rmsf_dataB = pt.rmsf(traj, mask=":58-107@CA byres")
rmsf_dataC = pt.rmsf(traj, mask=":110-162@CA byres")
rmsf_dataD = pt.rmsf(traj, mask=":165-219@CA byres")

residues=max(rmsf_dataA.T[0])
max_rmsf=max(rmsf_dataA.T[1])
np.savetxt('PLOTS/RMSF.csv', rmsf_data, fmt='%1.1f', delimiter='\t')
#Plot RMSF
plt.figure(2)                                                             

range_A= range(120,173)                                         # residue ranges of protein chains
range_B= range(124,172)
range_C= range(121,172)
range_D= range(120,173)

plt.plot(range_A, rmsf_dataA.T[1], label="Chain A",linewidth=0.7)
plt.plot(range_B, rmsf_dataB.T[1], label="Chain B",linewidth=0.7)
plt.plot(range_C, rmsf_dataC.T[1], label="Chain C",linewidth=0.7)
plt.plot(range_D, rmsf_dataD.T[1], label="Chain D",linewidth=0.7)


plt.xlabel('Residue', fontweight ='bold', fontsize=10)
plt.ylabel('RMSF (Angstrom)', fontweight ='bold', fontsize=10)
plt.xlim(119,173)
plt.ylim(0,max_rmsf+2)

plt.xticks((np.arange(120,173,4)),rotation=50, fontsize=8, labels=res[::4]) 
plt.yticks(fontsize=10)
plt.legend()




#_label_list = [1000,2000,3000,4000,5000]
#ax.set_xticks([0,200,400,600,800,1000])
#ax.set_xticklabels(x_label_list)

plt.savefig('PLOTS/RMSF_b.png', dpi=300, bbox_inches='tight', pad_inches=0.2)




