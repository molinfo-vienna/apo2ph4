import sys

from collections import Counter

import CDPL.Pharm as Pharm
import CDPL.Chem as Chem
import CDPL.MolProp as MolProp
        
ligand = Chem.BasicMolecule()
sdf_reader = Chem.MoleculeReader(sys.argv[1])

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
