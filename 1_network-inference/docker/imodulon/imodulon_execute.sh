#!/bin/bash

# if [ $# -le 1 ] || [ ! -f "${@:$#}" ]
# then
# 	../modulome-nextflow/4_optICA/run_ica.sh -h
# 	echo "FILE must come last."
# 	exit 1
# fi

# indir="$(dirname "${@:$#}")"

# datamash transpose < "${@:$#}" > "${indir}/expression_file.tsv"

# cd ../modulome-nextflow/4_optICA

# ./run_ica.sh "${@:1:$#}" "${indir}/expression_file.tsv"

# cd -

# rm "${indir}/expression_file.tsv"

python3 imodulon_execute2.py $@
