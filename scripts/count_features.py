import sys
import os.path as path
import CDPL.Pharm as Pharm
import CDPL.Math as Math
import CDPL.Chem as Chem
import CDPL.Base as Base
import CDPL.MolProp as MolProp
import numpy
from collections import Counter


def count_features():
    if len(sys.argv) < 2:
        print >> sys.stderr, 'Usage:', sys.argv[0], '[input.sdf]'
        sys.exit(2)
    ifs = Base.FileIOStream(sys.argv[1], 'r')

    ligand = Chem.BasicMolecule()
    sdf_reader = Chem.SDFMoleculeReader(ifs)

    lig_pharm = Pharm.BasicPharmacophore()
    pharm_gen = Pharm.DefaultPharmacophoreGenerator()
    ftr_list = list()

    while sdf_reader.read(ligand):
        Chem.perceiveSSSR(ligand, True)
        Chem.setAromaticityFlags(ligand, False)
        Chem.setRingFlags (ligand, False)
        MolProp.calcAtomHydrophobicities(ligand, False)
        
        pharm_gen.generate(ligand, lig_pharm)
        ftr_list += [Pharm.getType(ftr) for ftr in lig_pharm]
    
    with open("feature_count.txt", "w+") as writer:
        writer.write(str(dict(Counter(ftr_list))))
    
        
    


 


if __name__ == '__main__':
    count_features()
