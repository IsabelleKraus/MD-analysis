
# generate CHLD_cluster.. with MixMD Plugin (load grid pdb file) 

select pocket, 7r1o within 3 of CHLD_cluster_pts-grid
list=[]
iterate (pocket), list.append((resi, resn, name))
print list


