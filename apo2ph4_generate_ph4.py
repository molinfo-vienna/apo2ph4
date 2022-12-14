import sys
import os.path as path
import CDPL.Pharm as Pharm
import CDPL.Math as Math
import CDPL.Chem as Chem
import CDPL.Base as Base
import math
import argparse
from decimal import *
getcontext().prec = 6
from ast import literal_eval
import time
import numpy as np
from sklearn.cluster import DBSCAN
from copy import deepcopy


parser = argparse.ArgumentParser(description="Arguments for apo2ph4 pharmacophore generation")
req = parser.add_argument_group('required arguments')
req.add_argument('-i','--input_pml', type=str, required=True, help='path to input .pml file')
req.add_argument('-o','--output_pml', type=str, required=True, help='path to output .pml file')
req.add_argument('-g','--grid_maps', type=str, required=True, help='directory containing grid maps and feature_count.txt')
opt = parser.add_argument_group('optional arguments')
opt.add_argument('-n','--num_features', type=int, metavar='[int]', required=False, default=0, help='Number of desired total features, default=0 (auto)')
opt.add_argument('-p','--name', type=str,  metavar='\"Pharmacophore name\"', required=False, default='default', help='name of the pharmacophore model displayed in LigandScout')
opt.add_argument('-t','--distance_threshold', type=int, metavar='[int]', required=False, default=3, help='Distance threshold between features of the same type, default=3')
opt.add_argument('-v','--verbose', action='store_true', required=False, help='Verbose output (default: off)')
feat = parser.add_argument_group('features thresholds')
feat.add_argument('-H','--max_hydrophobic', type=int, metavar='[int]', required=False, default=4, help='Maximum number of hydrophobic features generated, default=4')
feat.add_argument('-D','--max_HBD', type=int, metavar='[int]', required=False, default=2, help='Maximum number of HBD features generated, default=2')
feat.add_argument('-P','--max_PI', type=int, metavar='[int]', required=False, default=1, help='Maximum number of PI features generated, default=1')
feat.add_argument('-N','--max_NI', type=int, metavar='[int]', required=False, default=1, help='Maximum number of NI features generated, default=1')
energy_threshold = parser.add_argument_group('energy thresholds')
energy_threshold.add_argument('--PI_energy', type=float, metavar='[float]', required=False, default=-0.7, help='Grid energy cutoff for PI features, default=-0.7 kcal/mol')
energy_threshold.add_argument('--NI_energy', type=float, metavar='[float]', default=0.9, help='Grid energy cutoff for NI features, default=0.9 kcal/mol')
energy_threshold.add_argument('--H_energy', type=float, metavar='[float]', default=-0.45, help='Grid energy cutoff for Hydrophobic features, default=-0.45 kcal/mol')
flags = parser.add_argument_group('flags for vector/planar features')
flags.add_argument('--ONE_FEATURE_PER_ORIGIN', type=int, choices=[0,1],required=False, default=1, help='Only allow one feature to be generated for each vector origin point, default=1')
flags.add_argument('--SUMMARIZE_OVERLAPPING_VECTORS', type=int, choices=[0,1],required=False, default=1, help='Summarizes features into a sphere if vectoral features from multiple origins are present at the same location, default=1')
flags.add_argument('--HBD_ONLY_SPHERES', type=int, choices=[0,1],required=False, default=0, help='Always display HBD as spheres, default=0, requires SUMMARIZE_OVERLAPPING_VECTORS=1')
flags.add_argument('--HBA_ONLY_SPHERES', type=int, choices=[0,1],required=False, default=0, help='Always display HBA as spheres, default=0, requires SUMMARIZE_OVERLAPPING_VECTORS=1')
flags.add_argument('--REMOVE_UNREALISIC_AR_ANGLES', type=int, choices=[0,1],required=False, default=1, help='Checks whether angle two adjacent aromatic features is realistic to be fulfilled by a real molecule, default=1')


args = parser.parse_args()

gridfile_path = args.grid_maps
if gridfile_path == "/": 
    gridfile_path = gridfile_path[:-1]


def process():
    start = time.time()
  
    ONLY_ONE_FEATURE_PER_ORIGIN = True #args.ONE_FEATURE_PER_ORIGIN   #should be True as it is unlikely a ligand satisfies a single HBD with more than one feature
    if args.ONE_FEATURE_PER_ORIGIN == 0: ONLY_ONE_FEATURE_PER_ORIGIN = False
    SUMMARIZE_OVERLAPPING_VECTORS = True   #should be True as important features could be missed due to individual scores of features
    if args.SUMMARIZE_OVERLAPPING_VECTORS == 0: SUMMARIZE_OVERLAPPING_VECTORS = False

    HBD_ONLY_SPHERES = False      #only valid if SUMMARIZE_OVERLAPPING_VECTORS is True
    if args.HBD_ONLY_SPHERES == 1: HBD_ONLY_SPHERES = True
    HBA_ONLY_SPHERES = False      #only valid if SUMMARIZE_OVERLAPPING_VECTORS is True
    if args.HBA_ONLY_SPHERES == 1: HBA_ONLY_SPHERES = True
    REMOVE_UNREALISIC_AR_ANGLES, MAX_AR_DISTANCE, AR_ANGLE_TOLERANCE = True, 4.0, 20 # Bool, Angstrom, Degrees;  filters out closeby aromatic features if their relative angle is unrealistic.
    if args.REMOVE_UNREALISIC_AR_ANGLES == 0: REMOVE_UNREALISIC_AR_ANGLES = False                                                                               
    
    
    try:    
        with open(gridfile_path + "/feature_count.txt", 'r') as feature_count:
           feature_count_dict = literal_eval(feature_count.readline())
    except FileNotFoundError: 
        print("feature_count.txt not found, run count_features.py first")
        sys.exit(2)
        

    def read_H_grid():
        newmap = open(args.grid_maps + '/receptor.C.map', 'r')
        lines = newmap.readlines()
        newmap.close()
        return lines

    def read_electrostatic_grid():
        newmap = open(args.grid_maps + '/receptor.e.map', 'r')
        lines = newmap.readlines()
        newmap.close()
        return lines

    def read_HBD_grid():
        newmap = open(args.grid_maps + '/receptor.HD.map', 'r')
        lines = newmap.readlines()
        newmap.close()
        return lines         

    def read_HBA_grid():
        newmap = open(args.grid_maps + '/receptor.OA.map', 'r')
        lines = newmap.readlines()
        newmap.close()
        return lines     

    def read_Ar_grid():
        newmap = open(args.grid_maps + '/receptor.A.map', 'r')
        lines = newmap.readlines()
        newmap.close()
        return lines    
        
    
    def get_average_grid_energy(gridfile, x_coord, y_coord, z_coord):
        spacing = Decimal(gridfile[3].split()[1])
        n_elements = list(map(Decimal,list(gridfile[4].split())[1:])) 
        center = list(map(Decimal,list(gridfile[5].split())[1:])) 
        origin = [center[0] - spacing * (n_elements[0]/2), center[1] - spacing * (n_elements[1]/2), center[2] - spacing * (n_elements[2]/2)]
        x_point = round((x_coord-float(origin[0]))/float(spacing))
        y_point = round((y_coord-float(origin[1]))/float(spacing))
        z_point = round((z_coord-float(origin[2]))/float(spacing))
        return float(gridfile[(z_point*(int(n_elements[1])+1)*(int(n_elements[0])+1)+y_point*(int(n_elements[0])+1)+x_point)+6])

    
    def count_psp_events(x_coord, y_coord, z_coord):
        spacing = Decimal(H_gridfile[3].split()[1])
        n_elements = list(map(Decimal,list(H_gridfile[4].split())[1:])) 
        center = list(map(Decimal,list(H_gridfile[5].split())[1:])) 
        origin = [center[0] - spacing * (n_elements[0]/2), center[1] - spacing * (n_elements[1]/2), center[2] - spacing * (n_elements[2]/2)]
        x_point = round((x_coord-float(origin[0]))/float(spacing))
        y_point = round((y_coord-float(origin[1]))/float(spacing))
        z_point = round((z_coord-float(origin[2]))/float(spacing))
        

        def count_x():
            positive_direction = False
            negative_direction = False
            x_min, x_max = 0, int(n_elements[0])
            x_current = x_point + 1
            for i in range(x_point, x_max + 1):
                if float(H_gridfile[(z_point*(int(n_elements[1])+1)*(int(n_elements[0])+1)+y_point*(int(n_elements[0])+1)+i)+6]) > 10:
                    positive_direction = True
                    break

            x_current = x_point - 1
            for i in range(x_point, x_min + 1, -1):
                if float(H_gridfile[(z_point*(int(n_elements[1])+1)*(int(n_elements[0])+1)+y_point*(int(n_elements[0])+1)+i)+6]) > 10:
                    negative_direction = True
                    break
            if positive_direction and negative_direction:
                return 1
            return 0


        def count_y():
            positive_direction = False
            negative_direction = False
            y_min, y_max = 0, int(n_elements[1])
            y_current = y_point + 1
            for i in range(y_point, y_max + 1):
                if float(H_gridfile[(z_point*(int(n_elements[1])+1)*(int(n_elements[0])+1)+i*(int(n_elements[0])+1)+x_point)+6]) > 10:
                    positive_direction = True
                    break

            y_current = y_point - 1
            for i in range(y_point, y_min + 1, -1):
                if float(H_gridfile[(z_point*(int(n_elements[1])+1)*(int(n_elements[0])+1)+i*(int(n_elements[0])+1)+x_point)+6]) > 10:
                    negative_direction = True
                    break
            if positive_direction and negative_direction:
                return 1
            return 0
                

        def count_z():      
            positive_direction = False
            negative_direction = False
            z_min, z_max = 0, int(n_elements[2])
            z_current = z_point + 1
            for i in range(z_point, z_max + 1):
                if float(H_gridfile[(i*(int(n_elements[1])+1)*(int(n_elements[0])+1)+y_point*(int(n_elements[0])+1)+x_point)+6]) > 10:
                    positive_direction = True
                    break

            z_current = z_point - 1
            for i in range(z_point, z_min + 1, -1):
                if float(H_gridfile[(i*(int(n_elements[1])+1)*(int(n_elements[0])+1)+y_point*(int(n_elements[0])+1)+x_point)+6]) > 10:
                    negative_direction = True
                    break
            if positive_direction and negative_direction:
                return 1
            return 0


        def count_diag_1():   # +x +y +z
            positive_direction = False
            negative_direction = False
            max_steps_pos = min(int(n_elements[0]) - x_point, int(n_elements[1]) - y_point, int(n_elements[2]) - z_point)
            max_steps_neg = -min(x_point, y_point, z_point)
            

            for i in range(0, max_steps_pos + 1):
                if float(H_gridfile[((z_point +i)*(int(n_elements[1])+1)*(int(n_elements[0])+1)+(y_point + i)*(int(n_elements[0])+1)+(x_point + i))+6]) > 10:       
                    positive_direction = True
                    break

            for i in range(0, max_steps_neg - 1, -1):
                if float(H_gridfile[((z_point  +i)*(int(n_elements[1])+1)*(int(n_elements[0])+1)+(y_point + i)*(int(n_elements[0])+1)+(x_point + i))+6]) > 10:
                    negative_direction = True
                    break
            if positive_direction and negative_direction:
                return 1
            return 0


        def count_diag_2():   # -x +y +z
            positive_direction = False
            negative_direction = False
            max_steps_pos = min(x_point, int(n_elements[1]) - y_point, int(n_elements[2]) - z_point)
            max_steps_neg = -min(int(n_elements[0]) - x_point, y_point, z_point)

            for i in range(0, max_steps_pos + 1):
                if float(H_gridfile[((z_point + i)*(int(n_elements[1])+1)*(int(n_elements[0])+1)+(y_point + i)*(int(n_elements[0])+1)+(x_point - i))+6]) > 10:       
                    positive_direction = True
                    break

            for i in range(0, max_steps_neg - 1, -1):
                if float(H_gridfile[((z_point  + i)*(int(n_elements[1])+1)*(int(n_elements[0])+1)+(y_point + i)*(int(n_elements[0])+1)+(x_point - i))+6]) > 10:
                    negative_direction = True
                    break
            if positive_direction and negative_direction:
                return 1
            return 0


        def count_diag_3():   # +x -y +z
            positive_direction = False
            negative_direction = False
            max_steps_pos = min(int(n_elements[0]) - x_point, y_point, int(n_elements[2]) - z_point)
            max_steps_neg = -min(x_point, int(n_elements[1]) - y_point, z_point)
            

            for i in range(0, max_steps_pos + 1):
                if float(H_gridfile[((z_point +i)*(int(n_elements[1])+1)*(int(n_elements[0])+1)+(y_point - i)*(int(n_elements[0])+1)+(x_point + i))+6]) > 10:       
                    positive_direction = True
                    break

            for i in range(0, max_steps_neg - 1, -1):
                if float(H_gridfile[((z_point  +i)*(int(n_elements[1])+1)*(int(n_elements[0])+1)+(y_point - i)*(int(n_elements[0])+1)+(x_point + i))+6]) > 10:
                    negative_direction = True
                    break
            if positive_direction and negative_direction:
                return 1
            return 0
 

        def count_diag_4():   # +x +y -z
            positive_direction = False
            negative_direction = False
            max_steps_pos = min(int(n_elements[0]) - x_point, int(n_elements[1]) - y_point, z_point)
            max_steps_neg = -min(x_point, y_point, int(n_elements[2]) - z_point)
            

            for i in range(0, max_steps_pos + 1):
                if float(H_gridfile[((z_point - i)*(int(n_elements[1])+1)*(int(n_elements[0])+1)+(y_point + i)*(int(n_elements[0])+1)+(x_point + i))+6]) > 10:       
                    positive_direction = True
                    break

            for i in range(0, max_steps_neg - 1, -1):
                if float(H_gridfile[((z_point  - i)*(int(n_elements[1])+1)*(int(n_elements[0])+1)+(y_point + i)*(int(n_elements[0])+1)+(x_point + i))+6]) > 10:
                    negative_direction = True
                    break
            if positive_direction and negative_direction:
                return 1
            return 0 
        return count_x() + count_y() + count_z() + count_diag_1() + count_diag_2() + count_diag_3() + count_diag_4()



    print('- Generating Pharmacophore')
    H_gridfile = read_H_grid()
    E_gridfile = read_electrostatic_grid()
    HBD_gridfile = read_HBD_grid()
    HBA_gridfile = read_HBA_grid()
    Ar_gridfile = read_Ar_grid()


    ifs = Base.FileIOStream(args.input_pml, 'r')

    pharm = Pharm.BasicPharmacophore()
    pml_reader = Pharm.PMLPharmacophoreReader(ifs)

    print('- Processing PML-file:', args.input_pml, '...') 
    feat_H_list, feat_AR_list, feat_NI_list, feat_PI_list, feat_HBD_list, feat_HBA_list, feat_EXCLU_list = [], [], [], [], [], [], set()
    while pml_reader.read(pharm):   
        for ftr in pharm:
            if Pharm.getType(ftr) == 1:
                ftr_pos = str(Chem.get3DCoordinates(ftr))[4:-1].split(',')
                energy = get_average_grid_energy(H_gridfile ,float(ftr_pos[0]), float(ftr_pos[1]), float(ftr_pos[2]))
                type_and_pos = [Pharm.getType(ftr), float(ftr_pos[0]), float(ftr_pos[1]), float(ftr_pos[2]), Pharm.getTolerance(ftr), Pharm.getGeometry(ftr), energy]
                feat_H_list.append(type_and_pos)    
            elif Pharm.getType(ftr) == 2:
                ftr_pos = str(Chem.get3DCoordinates(ftr))[4:-1].split(',')
                ftr_geo = str(Pharm.getOrientation(ftr))[4:-1].split(',')
                energy = get_average_grid_energy(Ar_gridfile ,float(ftr_pos[0]), float(ftr_pos[1]), float(ftr_pos[2]))
                type_and_pos = [Pharm.getType(ftr), float(ftr_pos[0]), float(ftr_pos[1]), float(ftr_pos[2]), Pharm.getTolerance(ftr), Pharm.getGeometry(ftr), float(ftr_geo[0]), float(ftr_geo[1]), float(ftr_geo[2]), energy]
                feat_AR_list.append(type_and_pos)
            elif Pharm.getType(ftr) == 3:
                ftr_pos = str(Chem.get3DCoordinates(ftr))[4:-1].split(',')
                energy = get_average_grid_energy(E_gridfile,float(ftr_pos[0]), float(ftr_pos[1]), float(ftr_pos[2]))                 
                type_and_pos = [Pharm.getType(ftr), float(ftr_pos[0]), float(ftr_pos[1]), float(ftr_pos[2]), Pharm.getTolerance(ftr), Pharm.getGeometry(ftr), energy]
                feat_NI_list.append(type_and_pos)
            elif Pharm.getType(ftr) == 4:
                ftr_pos = str(Chem.get3DCoordinates(ftr))[4:-1].split(',')
                energy = get_average_grid_energy(E_gridfile, float(ftr_pos[0]), float(ftr_pos[1]), float(ftr_pos[2]))
                type_and_pos = [Pharm.getType(ftr), float(ftr_pos[0]), float(ftr_pos[1]), float(ftr_pos[2]), Pharm.getTolerance(ftr), Pharm.getGeometry(ftr), energy]
                feat_PI_list.append(type_and_pos)
            elif Pharm.getType(ftr) == 5:
                ftr_pos = str(Chem.get3DCoordinates(ftr))[4:-1].split(',')
                ftr_geo = str(Pharm.getOrientation(ftr))[4:-1].split(',')
                ftr_length = Pharm.getLength(ftr)
                energy = get_average_grid_energy(HBD_gridfile,float(ftr_pos[0]), float(ftr_pos[1]), float(ftr_pos[2]))
                origin = [float(ftr_pos[0]) + Pharm.getLength(ftr) * float(ftr_geo[0]), float(ftr_pos[1]) + Pharm.getLength(ftr) * float(ftr_geo[1]), float(ftr_pos[2]) + Pharm.getLength(ftr) * float(ftr_geo[2])]
                type_and_pos = [Pharm.getType(ftr), float(ftr_pos[0]), float(ftr_pos[1]), float(ftr_pos[2]), Pharm.getTolerance(ftr), Pharm.getGeometry(ftr), float(ftr_geo[0]), float(ftr_geo[1]), float(ftr_geo[2]), origin, Pharm.getLength(ftr), energy]
                if energy < 0:
                    feat_HBD_list.append(type_and_pos)
            elif Pharm.getType(ftr) == 6:
                ftr_pos = str(Chem.get3DCoordinates(ftr))[4:-1].split(',')
                ftr_geo = str(Pharm.getOrientation(ftr))[4:-1].split(',')
                ftr_length = Pharm.getLength(ftr)
                energy = get_average_grid_energy(HBA_gridfile,float(ftr_pos[0]), float(ftr_pos[1]), float(ftr_pos[2]))
                origin = [float(ftr_pos[0]) - Pharm.getLength(ftr) * float(ftr_geo[0]), float(ftr_pos[1]) - Pharm.getLength(ftr) * float(ftr_geo[1]), float(ftr_pos[2]) - Pharm.getLength(ftr) * float(ftr_geo[2])]
                type_and_pos = [Pharm.getType(ftr), float(ftr_pos[0]), float(ftr_pos[1]), float(ftr_pos[2]), Pharm.getTolerance(ftr), Pharm.getGeometry(ftr), float(ftr_geo[0]), float(ftr_geo[1]), float(ftr_geo[2]), origin, Pharm.getLength(ftr), energy]
                if energy < 0:
                    feat_HBA_list.append(type_and_pos)
            elif Pharm.getType(ftr) == 7:
                ftr_pos = str(Chem.get3DCoordinates(ftr))[4:-1].split(',')
                type_and_pos = (Pharm.getType(ftr), float(ftr_pos[0]), float(ftr_pos[1]), float(ftr_pos[2]), Pharm.getTolerance(ftr), Pharm.getGeometry(ftr))    # a set is used here to throw out duplicates for performance reasons
                feat_EXCLU_list.add(type_and_pos)                                                               # touple used as set needs hashable value
            else:
                print("Feature of unknown type found") 
    if args.verbose:
        print("Number of features found from source .pml:")
        print("H", len(feat_H_list)) 
        print("Ar", len(feat_AR_list))  
        print("PI", len(feat_PI_list))  
        print("NI", len(feat_NI_list))  
        print("HBD", len(feat_HBD_list))  
        print("HBA", len(feat_HBA_list))  
        print("EXCL", len(feat_EXCLU_list))  
            

    feat_EXCLU_list = list(feat_EXCLU_list)                                                                     # convert set to list
    feat_EXCLU_list = [list(i) for i in feat_EXCLU_list]    
    threshold = args.distance_threshold # cutoff for scoring, everything beyond threshold contributes nothing to score)


    def score_and_filter(feature_list, threshold_multiplier = 1):        
        count = 0
        for i in feature_list:
            score = 0
            for j in feature_list:
                dist = math.sqrt(sum([(x1 - x2)**2 for (x1, x2) in zip(i[1:4], j[1:4])]))
                if dist < threshold: score += (1 + math.exp(4*(-1.0+dist)))** -1


            if feature_list[0][0] < 7:
                score = 100* (score)/(feature_count_dict.get(feature_list[0][0]))#100*(score - threshold)/(feature_count_dict.get(feature_list[0][0]))

            feature_list[count].insert(4, score) 
            count += 1



        feature_list.sort(key=lambda x: x[4], reverse=True)
        prox_matrix = [[math.sqrt(sum([(x1 - x2)**2 for (x1, x2) in zip(i[1:4], j[1:4])])) for j in feature_list] for i in feature_list]  # recalculating proximity matrix after sorting...
        filtered_list = []
        elements_in_list = []
        for i, j in enumerate(feature_list):
            if feature_list[0][0] == 4:
                adjacent_features = [k for k,l in enumerate(prox_matrix[i]) if l < threshold and feature_list[i][7] < PI_energy_cutoff]
            elif feature_list[0][0] == 3:
                adjacent_features = [k for k,l in enumerate(prox_matrix[i]) if l < threshold and feature_list[i][7] > NI_energy_cutoff]
            elif feature_list[0][0] == 1:
                adjacent_features = [k for k,l in enumerate(prox_matrix[i]) if l < threshold and feature_list[i][7] < H_energy_cutoff]
            else:
                adjacent_features = [k for k,l in enumerate(prox_matrix[i]) if l < threshold] 
            for m in adjacent_features:
                if m == i:   
                    filtered_list.append(j)
                    elements_in_list.append(i)
                    break
                if m in elements_in_list: break
                    
        return filtered_list




    def process_vector_features(feature_list):
        try:
            vector_array = np.empty((0,3), float)
            for i in feature_list:
                vector_array = np.vstack((vector_array, np.array(i[9])))    
           
            clustering = DBSCAN(eps=0.01, min_samples=3).fit(vector_array)
            arr_list = [[None] for _ in range(0, np.max(clustering.labels_)+1)]   
            [arr_list[j].append(i)  for i, j in zip(feature_list, clustering.labels_) if j >= 0]

            for i in arr_list: i.pop(0)

            unique, counts = np.unique(clustering.labels_, return_counts=True)
               
            feature_list = [i for j in [score_and_filter(i) for i in arr_list] for i in j]   
            feature_list.sort(key=lambda x: x[4], reverse=True)
        except: return feature_list

        if ONLY_ONE_FEATURE_PER_ORIGIN:
            origin_matrix = [[math.sqrt(sum([(x1 - x2)**2 for (x1, x2) in zip(i[10], j[10])])) for j in feature_list] for i in feature_list]       #proximity matrix of origins crds based on euclidean distance
            feature_list_temp = []
            elements_in_list = []
            for i, j in enumerate(feature_list):     #removes feature if feature of same origin was present in list beforehand)
                adjacent_features  = [k for k,l in enumerate(origin_matrix[i]) if l < 0.2]
                for m in adjacent_features:
                    if m == i:
                        feature_list_temp.append(j)
                        elements_in_list.append(i)
                        break
                    if m in elements_in_list: break
            feature_list = feature_list_temp 
                    

        if SUMMARIZE_OVERLAPPING_VECTORS:
            proximity_matrix = [[math.sqrt(sum([(x1 - x2)**2 for (x1, x2) in zip(i[1:4], j[1:4])])) for j in feature_list] for i in feature_list]  #proximity matrix of ftr crds based on euclidean distance
            feature_list_temp = []
            elements_in_list = []
            for i, j in enumerate(feature_list):
                adjacent_features = [k for k,l in enumerate(proximity_matrix[i]) if l < 1.5]
                for m in adjacent_features:
                    if m == i:
                        j[4] = sum([feature_list[m][4] for m in adjacent_features])
                        if len(adjacent_features) > 1: j[6] = 1
                        feature_list_temp.append(j)
                        elements_in_list.append(i)
                        break
                    if m in elements_in_list: break
            feature_list = feature_list_temp

        return feature_list

    def process_aromatic_features(feature_list):
            start2 = time.time()
            vector_array = np.empty((0,3), float)#

            for i in feature_list:
                vector_array = np.vstack((vector_array, np.array(i[6:9]))) 

            try:   
                clustering = DBSCAN(eps=0.15, min_samples=2).fit(vector_array)
                arr_list = [[None] for _ in range(0, np.max(clustering.labels_)+1)]   
                [arr_list[j].append(i)  for i, j in zip(feature_list, clustering.labels_) if j >= 0]
            except: return feature_list

            for i in arr_list: i.pop(0)

            unique, counts = np.unique(clustering.labels_, return_counts=True)

            feature_list = [i for j in [score_and_filter(i) for i in arr_list] for i in j]   
            feature_list.sort(key=lambda x: x[4], reverse=True)

            proximity_matrix = [[math.sqrt(sum([(x1 - x2)**2 for (x1, x2) in zip(i[1:4], j[1:4])])) for j in feature_list] for i in feature_list]  #proximity matrix of ftr crds based on euclidean distance
            feature_list_temp = []
            elements_in_list = []
            for i, j in enumerate(feature_list):
                adjacent_features = [k for k,l in enumerate(proximity_matrix[i]) if l < 1.5]
                for m in adjacent_features:
                    if m == i:
                        j[4] = sum([feature_list[m][4] for m in adjacent_features])             #while it seems like, that this defeats the whole purpose of clustering, this at least insures, that
                        if len(adjacent_features) > 1: j[6] = 1                                 #a highly populated plane angle is selected  (an outlier could be selected by chance, it also removes outliers
                        feature_list_temp.append(j)                                             #that do not have a neighbor within ~10° and thus also affects the scoring
                        elements_in_list.append(i)
                        break
                    if m in elements_in_list: break
            feature_list = feature_list_temp

            if REMOVE_UNREALISIC_AR_ANGLES:
                angle_tolerance_distance = math.sqrt(2-2*math.cos(math.radians(AR_ANGLE_TOLERANCE)))
                filtered_out = []
                for l, i in enumerate(feature_list):
                    for k, j in enumerate(feature_list):
                        if math.sqrt(sum([(x1 - x2)**2 for (x1, x2) in zip(i[1:4], j[1:4])])) <= MAX_AR_DISTANCE and math.sqrt(sum([(x1 - x2)**2 for (x1, x2) in zip(i[6:9], j[6:9])])) > angle_tolerance_distance and k > l and k not in filtered_out:
                            filtered_out.append(j)
                            if args.verbose: print("Filtered out Aromatic feature", j, "as it did exceed the defined angle tolerance next to other aromatic feature")
                feature_list = [i for i in feature_list if i not in filtered_out]
                            
            return feature_list


    feature_dict = {
        "1" : "H  ",
        "2" : "AR ",
        "3" : "NI ",
        "4" : "PI ",
        "5" : "HBD",
        "6" : "HBA",
        "7" : "EXC"
    }

    
    exclusion_volumes = score_and_filter(feat_EXCLU_list)
  

    PI_energy_cutoff = args.PI_energy
    NI_energy_cutoff = args.NI_energy
    H_energy_cutoff = args.H_energy
    curr_NI, max_NI = 0, args.max_NI 
    curr_PI, max_PI = 0, args.max_PI
    curr_H, max_H = 0, args.max_hydrophobic
    curr_HBD, max_HBD = 0, args.max_HBD 


    all_features_filtered = score_and_filter(feat_H_list) + process_aromatic_features(feat_AR_list) + process_vector_features(feat_HBA_list) + process_vector_features(feat_HBD_list) + score_and_filter(feat_NI_list) + score_and_filter(feat_PI_list) 
 


    for j, i in enumerate(all_features_filtered):     
        if i[0] == 5: all_features_filtered[j][4] = i[4] * 0.6  

    all_features_filtered.sort(key=lambda x: x[4], reverse=True)
    


    pharmacophore_list = []
    

    for i in all_features_filtered:
        if i[0] == 1:
            if i[-1] < curr_H < max_H:
                pharmacophore_list.append(i)
                curr_H += 1
        elif i[0] == 3:
            if i[-1] >  curr_NI < max_NI:
                pharmacophore_list.append(i)
                curr_NI += 1
        elif i[0] == 4:
            if i[-1] < curr_PI < max_PI:
                pharmacophore_list.append(i)
                curr_PI += 1
        elif i[0] == 5:
            if curr_HBD < max_HBD:
                pharmacophore_list.append(i) 
                curr_HBD += 1         
        elif i[0] == 6:
            pharmacophore_list.append(i)   
        else: pharmacophore_list.append(i)


    score_top_four_avg = 0
    for i in range(0,4):
        score_top_four_avg += pharmacophore_list[i][4]/4

    max_pharmacophores = 0    
    default_max_pharmacophores = args.num_features
    if default_max_pharmacophores != 0:
        max_pharmacophores = default_max_pharmacophores
    else:
        try:  #for the rare case of the last pharmacophore feature still fitting the threshold
            while pharmacophore_list[default_max_pharmacophores][4] > 0.2*score_top_four_avg:
                default_max_pharmacophores += 1  
            max_pharmacophores = default_max_pharmacophores
        except: max_pharmacophores = default_max_pharmacophores



    

    print("- Done!")
    print("\nNumber of features generated: ", max_pharmacophores)
    print("Generated pharmacophore features:")
    print("Type\t  X\t  Y\t  Z\tScore\tBuriedness(PSP)\t\tGrid energy [kcal/mol]")
    try:    
        for i in range(0,max_pharmacophores):
            print(feature_dict[str(pharmacophore_list[i][0])], "\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t".format(pharmacophore_list[i][1], pharmacophore_list[i][2], pharmacophore_list[i][3], pharmacophore_list[i][4]), count_psp_events(pharmacophore_list[i][1], pharmacophore_list[i][2], pharmacophore_list[i][3]), "\t\t\t{:.3}".format(pharmacophore_list[i][-1]))
    except: pass



    def write_pharm():
        new_pharm = Pharm.BasicPharmacophore()
        feat_set = Pharm.FeatureSet()           
        for i in pharmacophore_list[:max_pharmacophores]: 
            new_feature = new_pharm.addFeature()
            Pharm.setType(new_feature,i[0])
            coords = Math.Vector3D()
            coords = (i[1], i[2], i[3])
            Chem.set3DCoordinates(new_feature,coords)
            if i[0] == 2:
                Pharm.setTolerance(new_feature, 1.2)    
                geo_coords = Math.Vector3D()
                geo_coords = (i[7], i[8], i[9])
                Pharm.setOrientation(new_feature,geo_coords)
                Pharm.setGeometry(new_feature, 3)
            elif i[0] == 5 or i[0] == 6:
                geo_coords = Math.Vector3D()
                geo_coords = (i[7], i[8], i[9])
                Pharm.setOrientation(new_feature,geo_coords)
                Pharm.setTolerance(new_feature, i[5]) 
                Pharm.setLength(new_feature, i[-2])
                Pharm.setGeometry(new_feature, i[6])
                if HBD_ONLY_SPHERES and SUMMARIZE_OVERLAPPING_VECTORS and i[0] == 5:
                    Pharm.setGeometry(new_feature, 1)
                if HBA_ONLY_SPHERES and SUMMARIZE_OVERLAPPING_VECTORS and i[0] == 6:
                    Pharm.setGeometry(new_feature, 1)
            else: 
                Pharm.setTolerance(new_feature, i[5]*1) 
                Pharm.setGeometry(new_feature,1)
            feat_set.addFeature(new_feature)

        for i in exclusion_volumes:
            new_feature = new_pharm.addFeature()
            Pharm.setType(new_feature,i[0])
            coords = Math.Vector3D()
            coords = (i[1], i[2], i[3])
            Pharm.setTolerance(new_feature, 1) 
            Chem.set3DCoordinates(new_feature,coords)
            feat_set.addFeature(new_feature)
        
        if args.name == "default":
            Pharm.setName(feat_set, ''.join((str(max_pharmacophores), 'ftr')))
        else:
            Pharm.setName(feat_set, args.name)
        writer = Pharm.PMLFeatureContainerWriter(Base.FileIOStream(args.output_pml,'w'))
        writer.write(feat_set)
    
    write_pharm()
    print ("\ntime elapsed:", round(time.time()-start,1), "seconds")


if __name__ == '__main__':
    process()
