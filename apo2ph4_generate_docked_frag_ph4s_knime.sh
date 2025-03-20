#!/bin/bash

SCRIPT=$(readlink -f $0)
SCRIPTPATH=`dirname $SCRIPT`

mkdir -p pharmacophores

knime --launcher.suppressErrors -nosave -reset -nosplash -application org.knime.product.KNIME_BATCH_APPLICATION -workflowDir="${SCRIPTPATH}/knime_fragment_SB_pharmacophore_gen/"  -workflow.variable=var_output,"./pharmacophores/pharmacophores.pml",String -workflow.variable=var_input,"./complex_list.list",String -vmargs -Xmx32000m -Xms10048m 
