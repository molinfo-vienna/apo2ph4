import sys



mw = {
    'H' : 1,
    'B' : 10.8,
    'C' : 12,
    'N' : 14,
    'O' : 16,
    'F' : 19,
    'SI': 28.1,
    'Si': 28.1,
    'P' : 31,
    'S' : 32.1,
    'CL': 35.5,
    'Cl': 35.5,
    'BR': 79.9,
    'Br': 79.9,
    'I' : 126.9,
    'PS': 1
}



with open(sys.argv[1], 'r') as f:
    a = []
    line = f.readline()
    total_weight = 0
    x_sum, y_sum, z_sum = 0, 0, 0
    while line:
        if line.startswith('HETATM'):
            weight = mw[line[76:78].strip()]   
            total_weight += weight
            x_sum += float(line[30:38].strip())*weight
            y_sum += float(line[38:46].strip())*weight
            z_sum += float(line[46:54].strip())*weight
        line = f.readline()
    x = round(x_sum/total_weight, 2)
    y = round(y_sum/total_weight, 2)
    z = round(z_sum/total_weight, 2)


with open('vinaconf.txt', 'w+') as conf:
    conf.write('receptor = receptor.pdbqt\n')
    conf.write('center_x = ')
    conf.write(str(x))
    conf.write('\ncenter_y = ')
    conf.write(str(y))
    conf.write('\ncenter_z = ')
    conf.write(str(z))
    conf.write('\nsize_x = 20')
    conf.write('\nsize_y = 20')
    conf.write('\nsize_z = 20')
    conf.write('\nexhaustiveness = 20')
    conf.write('\ncpu = 8')
    conf.write('\nnum_modes = 2')
    conf.write('\nenergy_range = 2')


