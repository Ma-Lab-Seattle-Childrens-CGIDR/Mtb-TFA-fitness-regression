#!/bin/bash

# go straight to help text if there aren't at least 2 args and the last arg isn't a real file
if [ $# -le 1 ] || [ ! -f "${@:$#}" ]
then
	cmonkey2 -h
	echo "ratios file must come last."
	exit 1
fi

indir="$(dirname "${@:$#}")"
run="$(date +%s)"

datamash transpose < "${@:$#}" > "${indir}/expression_file-${run}.tsv"

# pass all received arguments minus the last, and add this temporary file at the end
# this is because cmonkey2 expects a spreadsheet in transposed format compared to the other methods
cmonkey2 "${@:1:$#}" "${indir}/expression_file-${run}.tsv"

rm "${indir}/expression_file-${run}.tsv"
