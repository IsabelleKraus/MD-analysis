parm z
trajin FtuFabI.nc
reference FtuFabI.WT.pdb
nativecontacts :7,8,9,11,14,21,29,31 :58-61 \   # [('7', 'ASN'), ('8', 'VAL'), ('9', 'ILE'), ('11', 'ASP'), ('14', 'ASN'), ('21', 'ARG'), ('29', 'ASP'), ('31', 'ASP')]
                                                # resi 58-61 → Ligand 

byresidue out nc.all.res.dat mindist maxdist \  
distance 3.0 reference map mapout resmap.gnu \
contactpdb Loop-NDP.pdb \
series seriesout native.dat

 .%OS .%OH distance 8.0 contact.out byresidue map mapout map.out 
