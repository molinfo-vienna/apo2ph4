##
# gen_ia_ph4s.py 
#
# Copyright (C) 2025 Thomas Seidel <thomas.seidel@univie.ac.at>
#
##

import sys

import CDPL.Chem as Chem
import CDPL.Biomol as Biomol
import CDPL.Pharm as Pharm
import CDPL.MolProp as MolProp


# reads and preprocesses the specified docked fragment complex
def readAndPrepareComplex(compl: str, compl_mol: Chem.Molecule) -> None:
    reader = Chem.MoleculeReader(compl) # create reader for complex (format specified by file extension)
    
    try:
        if not reader.read(compl_mol):  # read complex
            sys.exit(f'Error: reading complex file \'{compl}\' failed')

    except Exception as e:
        sys.exit(f'Error: reading complex file \'{compl}\' failed: {str(e)}')

    # preprocess the read complex structure
    try:
        Chem.perceiveSSSR(compl_mol, False)
        Chem.setRingFlags(compl_mol, False)
        Chem.calcImplicitHydrogenCounts(compl_mol, False);
        Chem.perceiveHybridizationStates(compl_mol, False);
        Chem.setAromaticityFlags(compl_mol, False);

        if Chem.makeHydrogenComplete(compl_mol):                    # make implicit hydrogens (if any) explicit
            Chem.calcHydrogen3DCoordinates(compl_mol)               # calculate 3D coordinates for the added expl. hydrogens
            Biomol.setHydrogenResidueSequenceInfo(compl_mol, False) # set residue information for the added expl. hydrogens

        MolProp.calcAtomHydrophobicities(compl_mol, False)          # calculate atom hydrophobicity values (needed for hydrophobic
                                                                    # pharm. feature generation)
    except Exception as e:
        sys.exit(f'Error: processing of complex file \'{compl}\' failed: {str(e)}')

# split complex structure into separate receptor and ligand structures
def splitComplex(compl_mol: Chem.Molecule) -> tuple:
    rec_mg = Chem.Fragment()
    lig_mg = Chem.Fragment()
    first_atom = None
    
    for atom in compl_mol.atoms:
        if not first_atom:
            lig_mg.addAtom(atom)
            first_atom = atom
            continue

        if Biomol.getResidueCode(atom) == Biomol.getResidueCode(first_atom) and Biomol.getResidueSequenceNumber(atom) == Biomol.getResidueSequenceNumber(first_atom):
            lig_mg.addAtom(atom)
        else:
            rec_mg.addAtom(atom)

    for bond in compl_mol.bonds:
        if lig_mg.containsAtom(bond.getBegin()) and lig_mg.containsAtom(bond.getEnd()):
            lig_mg.addBond(bond)
        else:
            rec_mg.addBond(bond)

    Chem.extractSSSRSubset(compl_mol, lig_mg, True)
    Chem.extractSSSRSubset(compl_mol, rec_mg, True)
    
    return (rec_mg, lig_mg)
            
def main() -> None:
    compl_mol = Chem.BasicMolecule()                    # create an instance of the default implementation of the Chem.Molecule interface
    ph4_writer = Pharm.FeatureContainerWriter(sys.argv[2]) # create writer for the generated pharmacophores (format specified by file extension)
 
    ia_ph4 = Pharm.BasicPharmacophore()                 # create an instance of the default implementation of the Pharm.Pharmacophore interface
    ph4_gen = Pharm.InteractionPharmacophoreGenerator() # create an instance of the pharmacophore generator

    ph4_gen.addExclusionVolumes(True)                   # specify whether to generate exclusion volume spheres on pharm. feature atoms of interacting residues

    # disable unsupported feature types
    ph4_gen.envPharmacophoreGenerator.enableFeature(Pharm.FeatureType.HALOGEN_BOND_ACCEPTOR, False) 
    ph4_gen.corePharmacophoreGenerator.enableFeature(Pharm.FeatureType.HALOGEN_BOND_DONOR, False)
        
    for compl in open(sys.argv[1], 'r').readlines():
        compl = compl.strip()
        
        print(f'- Processing complex \'{compl}\'')
        
        readAndPrepareComplex(compl, compl_mol)   # read and preprocess the comples
        
        try:
            rec_mg, lig_mg = splitComplex(compl_mol) # split complex into receptor and ligand
            
            ph4_gen.generate(lig_mg, rec_mg, ia_ph4, True) # generate the pharmacophore (True = extract ligand environment residues on-the-fly)

            try:
                if not ph4_writer.write(ia_ph4): # output pharmacophore
                    sys.exit(f'Error: writing interaction pharmacophore of complex \'{compl}\' failed')

            except Exception as e:               # handle exception raised in case of severe write errors
                sys.exit(f'Error: writing interaction pharmacophore of complex \'{compl}\' failed: {str(e)}')
                
        except Exception as e:                   # handle exception raised in case of severe processing errors
            sys.exit(f'Error: interaction pharmacophore generation for complex \'{compl}\' failed: {str(e)}')

    print('- Done!')

    ph4_writer.close()
    sys.exit(0)
        
if __name__ == '__main__':
    main()
