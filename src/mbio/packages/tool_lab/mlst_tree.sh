python $1 $2 martix_result.xls
rapidnj martix_result.xls -i pd -o t > mlst_tree.nwk
sed "s/'//g" mlst_tree.nwk -i
