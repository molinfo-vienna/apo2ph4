from pymol import cmd
import os
import sys
pathname = os.path.dirname(sys.argv[0])   

cmd.feedback('disable', 'selector', 'error')


if len(sys.argv) < 5:
    print("Usage apo2ph4_definine_binding_site.py [receptor.PDB/CIF etc] [x-coordinate] [y-coordinate] [z-coordinate]")
    sys.exit(2)



cmd.load(sys.argv[1])
cmd.remove('organic')
cmd.remove('solvent')
cmd.remove('inorganic')
cmd.load(pathname+"/other/benzene.sdf")
expression = "(x,y,z)=(x+" + str(sys.argv[2]) + ",y+" + str(sys.argv[3]) + ",z+" + str(sys.argv[4]) + ")"
cmd.alter_state("1", "benzene", expression)
cmd.save(sys.argv[1][:-4]+"_prepared.pdb")


cmd.quit()
