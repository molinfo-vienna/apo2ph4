# Instructions for the generation of *apo* site pharmacophores by the [Apo2ph4](https://pubs.acs.org/doi/10.1021/acs.jcim.2c00814) method

Copyright © 2022 Jörg Heider and Thomas Seidel

**Requirements**

- Python 3
- CDPKit python bindings (<https://pypi.org/project/CDPKit>)
- scikit-learn python package (<https://pypi.org/project/scikit-learn>)
- Pymol python bindings (e.g. <https://pymol.org/conda>)
- MGLTools (1.5.6 or newer) (<https://ccsb.scripps.edu/mgltools/downloads>)
- OpenBabel (<http://openbabel.org/docs/Installation/install.html>)
- KNIME (optional - see below, 4.1.2 or newer, <https://www.knime.com/downloads>)
- LigandScout KNIME exensions (optional - see below, available via KNIME update site https://www.inteligand.com/knime-nodes)
- AutoDock4 including AutoGrid4 (<https://autodock.scripps.edu/download-autodock4>)
- AutoDock Vina (<https://vina.scripps.edu/downloads>)

**Required to be accessible via the `PATH` environment variable:**

- MGLTools Utilities
- OpenBabel
- KNIME (optional - see below)
- AutoDock4
- AutoDock Vina

**Recommended way to install required external software:**

A quite simple and proven to work possibility to install all external software
packages (except *KNIME* whose use is optional) needed by the *Apo2ph4* workflow scripts
is to create a dedicated [conda](https://github.com/conda-forge/miniforge) environment and perform an installation of the required
packages by executing the commands below (assuming an active conda base environment)
in the given order:

```console
conda config --add channels conda-forge # already added when using Miniforge conda
conda config --add channels bioconda
conda create -n apo2ph4 python==3.10
conda activate apo2ph4
conda install pymol-open-source
conda install openbabel
conda install autodock
conda install autogrid
conda install autodock-vina
conda install mgltools
pip3 install cdpkit
pip3 install scikit-learn
```
The above environment setup commands need to be executed only once. The
*Apo2ph4* scripts will then run without any issues when executed from within an
active `apo2ph4` environment.

### Quickstart guide

The workflow comprises four steps (scripts) where one of them is
optional. If desired, these scripts may be executed in succession using
a bash script to require no intervention after a binding pocket is
selected.

**1\) apo2ph4\_define\_binding\_site.py** (Optional)

If a binding site is defined by specifying cartesian coordinates
`apo2ph4_define_binding_site.py` must be run beforehand. It takes
one PDB file and the x, y and z coordinates as positional arguments. This
script removes all ligands, ions, and solvent molecule before placing a
dummy molecule at the defined coordinates.

For example:

```console
python3 apo2ph4_define_binding_site.py [protein].pdb 141.021 55.154 89.114
```
This will generate *\[protein\]\_prepared.pdb* in the same folder as the
input PDB.

This file may also be opened in *LigandScout* which now recognizes the
binding site due to placement of a dummy ligand molecule.
 
**2\) apo2ph4\_prepare\_and\_dock.sh**

The second step of the workflow is a *bash* script that takes two
positional arguments:

```console
bash apo2ph4_prepare_and_dock.sh [fragment_library].sdf [protein]_prepared.pdb
```
*\[fragment\_library\].sdf* should contain the fragment database with
3D-coordinates in SDF format (this database may for example be prepared
by generating 3D-coordinates using the tool [structgen](https://cdpkit.org/applications/structgen.html)
and selecting SDF as output format or one can use the file
*fragments.sdf* located in the *data* directory).

*\[protein\]\_prepared.pdb* was either prepared beforehand (step 1) or
should be a PDB file only containing a single ligand that defines the
binding site (other ligands, ions and solvent molecules must have been
removed). The script should ideally be called from the directory that
contains the PDB file .

The necessary grid energy files are placed in the parent
folder as \*.map (e.g. *receptor.C.map*) files which also contains the
*feature\_count.txt* file required for subsequent steps.

*Vina* docking parameters are hardcoded, these may be changed by modifying
the `write_vina_conf.py` file in the workflow *scripts*
directory.

Intermediate files such as docked poses, mock ligand/protein complexes
are stored in directory *./tempdock*.

**Notes:**

The contents of *./tempdock* are not required anymore after the workflow step 3
has finished and are only kept for debugging purposes. As this data
files consume a lot of space (~500 MB for 200 fragments) it is recommended
to delete *./tempdock* (after step 3!) should these files not be desired.

**3a\) apo2ph4\_generate\_docked\_frag\_ph4s\_knime.sh**

This script generates interaction pharmacophores of the docked fragments
combined into a single *LigandScout* \*.pml file named
*./pharmacophores/pharmacophores.pml* in the current working
directory by means of a *KNIME* workflow executed in batch mode.

Script execution:

```console
bash apo2ph4_generate_docked_frag_ph4s_knime.sh
```
**Notes:**

*KNIME* will be executed with pre-set memory parameters.
If not enough resources are available or for optimal
performance the following parameters in the last line of
`apo2ph4_generate_docked_frag_ph4s_knime.sh` should be modified: *-vmargs
-Xmx32000m -Xms10048m*

*-Xmx32000m* defines the maximum available memory (in this case 32 GB) and
*-Xms10048m* defines the initial memory reserved (in this case 10 GB)

Pharmacophore generation is carried out using a subset of the KNIME nodes provided by the *LigandScout KNIME exensions*.
The *LigandScout* exensions are proprietary software and a valid license will be required to be able to successfully run
the KNIME workflow!

**3b\) apo2ph4\_generate\_docked\_frag\_ph4s\_cdpkit.sh**

This variant uses *CDPKit* functionality to generate the interaction pharmacophores of the docked fragments.
Since *CDPKit* is free software, no proprietary software license will be required to sucessfully run the script.

Script execution:

```console
bash apo2ph4_generate_docked_frag_ph4s_cdpkit.sh
```

The individual pharmacophores are as well written to a single file named
*./pharmacophores/pharmacophores.pml*. The *apo* pharmacophores
resulting from the last step of the overall workflow will not be 100% equal to the ones generated from
*LigandScout* fragment pharmacophores but are of comparable quality and will neither be worse nor better
on average.

**4\) apo2ph4\_generate\_ph4.py**

The last step of the workflow generates a final pharmacophore model and
requires the *pharmacophores.pml*, *grid map* and
*feature\_count.txt* files generated by the previous steps.

The number of pharmacophore features, energy thresholds as well as other
settings may be specified by optional parameters.

Execute `apo2ph4_generate_ph4.py -h` to see all available parameters
for advanced usage.

Basic example (mandatory arguments only):

```console
python3 apo2ph4_generate_ph4.py -i pharmacophores/pharmacophores.pml -o pharmacophores/my_ph4.pml -g ./
```
Advanced example:

```console
python3 apo2ph4_generate_ph4.py -i pharmacophores/pharmacophores.pml -o pharmacophores/my_ph4.pml -g ./ -n 8 --name "my pharmacophore model" -H 5 --H_energy -0.6 -v
```
It is recommended to examine the final pharmacophore model visually by
means of *LigandScout* to make sure no features have been placed outside
the actual binding pocket.

The resulting *apo* pharmacophores may then be used as query in pharmacophore-based
virtual screening campaigns (using [psdscreen](https://cdpkit.org/applications/psdscreen.html) or
*LigandScout*) or as source of valuable information for lead compound optimization.
