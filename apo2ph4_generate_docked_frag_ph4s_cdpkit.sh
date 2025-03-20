#!/bin/bash

SCRIPT=$(readlink -f $0)
SCRIPTPATH=`dirname $SCRIPT`

mkdir -p pharmacophores

python3 ${SCRIPTPATH}/scripts/gen_ia_ph4s.py ./complex_list.list ./pharmacophores/pharmacophores.pml
