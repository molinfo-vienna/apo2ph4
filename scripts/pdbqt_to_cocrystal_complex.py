from pymol import cmd
import os

cmd.feedback('disable', 'selector', 'error')

[cmd.load(pdbqt) for pdbqt in sorted(os.listdir('.')) if pdbqt.endswith("docked.pdbqt")]

cmd.split_states('*')

[cmd.remove(i) for i in cmd.get_object_list() if 'docked_' not in i]
[cmd.alter(j, ''.join(('resn=\'', str(i), '\''))) for i, j in enumerate(cmd.get_object_list())]
cmd.load('receptor.pdbqt')


[cmd.save(''.join(('fragment_complex', str(i), '.pdb')), ''.join(('receptor or ', j))) for i, j in enumerate(cmd.get_object_list()[:-1])]




with open ('../complex_list.list', 'w+') as complist: [complist.write(os.getcwd() + '/' + i + '\n') for i in sorted(os.listdir('.')) if i.startswith('fragment_complex')]
    
cmd.quit()
