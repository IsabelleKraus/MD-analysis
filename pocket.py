
# generate CHLD_cluster.. with MixMD Plugin (load grid pdb file) 

from pymol import cmd

def findBindingPocketResidues(selection="protein", cluster="CHLD_cluster", cutoff=3.0, doShow=0, quiet=1):
    """
DESCRIPTION
    Finds those residues of a possible binding pocket of a protein
    that are within 'cutoff' of an occupancy area found by MixMD.
USAGE
    findBindingPocketResidues [ selection, [cluster, [ cutoff, [ doShow ]]]]
ARGUMENTS
    selection = string: object or selection in which to find exposed
    residues {default: protein}
    cluster = string: object or selection in which to find high occupancy of cosolvents
    cutoff = float: cutoff of what is in binding pocket or not {default: 3.0 Ang**2}
RETURNS
    (list: (resi, resn ) )
        A Python list of residue numbers and names corresponding
        to those residues within the cutoff.
    """
    cutoff, doShow, quiet = float(cutoff), int(doShow), int(quiet)

    selName = findBindingPocketResidues(selection, cluster, cutoff, quiet)

    cmd.select(pocket, byres(selection) within cutoff of cluster)
    list=[]
    cmd.iterate(pocket and name CA), list.append((resi, resn))
    print list
    
    


