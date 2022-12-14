#arg1 = path to ligands sdf
#arg2 = path to receptor path to receptor with  core ligand
#!/bin/bash

SCRIPT=$(readlink -f $0)
SCRIPTPATH=`dirname $SCRIPT`

rm tempdock -rf
mkdir tempdock
mkdir pharmacophores
grep -E '(HETATM|CONECT)' $2 >> tempdock/coreligand.pdb   #extracts ligand into separate pdb
grep -vE '(HETATM|CONECT)' $2 >> tempdock/receptor.pdb     # extracts apo protein into separate pdb



babel $1 tempdock/ligand.pdb -m    #converts ligands to individual .pdb files

source ${SCRIPTPATH}/venv/bin/activate  # goes into python2 environment, needed for preparing pdbqt files
lignum=$(ls -l tempdock/ligand*.pdb | wc -l) 
cd tempdock
prepare_receptor4.py -r receptor.pdb -o receptor.pdbqt   #prepares receptor for autogrid and for docking
prepare_ligand4.py -l coreligand.pdb -o coreligand.pdbqt  #need this to generate autogrid file
prepare_gpf4.py -l coreligand.pdbqt -r receptor.pdbqt -y  -p  npts='100,100,100' -o out.gpf  -p ligand_types='C, OA, HD, A' -p spacing='0.300'  #prepares config file for autogrid4

for (( i = 1; i<=$lignum; i++ ))   
    do
        prepare_ligand4.py -l ligand${i}.pdb #this is not neccessary since we use the same ligands all the time, might remove...
        echo "processed ligand $i"
    done  #for loop that converts all ligand.pdb files in ligand.pdbqt files needed for docking

deactivate #leaves virtual environment



python3 ${SCRIPTPATH}/scripts/write_vina_conf.py ../$2  #writes vina configuration file
autogrid4 -p out.gpf  #path to autogrid executable, creates grid files needed later on
mv *.map ../




for (( i = 1; i<=$lignum; i++ ))
    do
    ligand=$(ls ligand${i}.pdbqt)
    echo "docking ligand $ligand"
    vina --config ./vinaconf.txt --ligand $ligand --out ligand${i}_docked.pdbqt
    done  #submits every ligand to docking



python3 ${SCRIPTPATH}/scripts/pdbqt_to_cocrystal_complex.py  #creates pseudo cocystal structures for kniome + a list to use as knime input
cd ..
python3 ${SCRIPTPATH}/scripts/count_features.py $1   #



knime --launcher.suppressErrors -nosave -reset -nosplash -application org.knime.product.KNIME_BATCH_APPLICATION -workflowDir="${SCRIPTPATH}/knime_fragment_SB_pharmacophore_gen/"  -workflow.variable=var_output,"./pharmacophores/pharmacophores.pml",String -workflow.variable=var_input,"./complex_list.list",String -vmargs -Xmx32000m -Xms10048m 
