#!/bin/bash

# SCRIPT PARA EXTRAER LAS TASAS DE MAPPING DE LOS ARCHIVOS .out DE STAR
cd refseq_k15
grep '\[ mates1 \]' *.out | sed -E 's|.*\/([^/]+)_trim_.*|\1|' > names.txt

# Extrae el mapping rate
grep 'Mapping rate' *.out | awk -F'= ' '{print $2}' | sed 's/ ->.*//' > rates.txt

# Combina ambos en un archivo final
echo -e "Samples\tMappingRate" > mapping_summary.tsv
paste names.txt rates.txt >> mapping_summary.tsv
