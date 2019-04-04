#/bin/bash

parameter=$1
## Backup jedi.params
cp jedi.params jedi.params.backup
# Call jedi-paramscan.py
python jedi-paramscan.py -p $parameter
# Perform the scan of the parameter
rm jedi.params

echo "Scaning JEDI using SITE keyword"
for params in jedi-scan-*.params
do
	foldername=`basename $params '.params'`
	mkdir $foldername
	ln -s $params jedi.params
	cp $params $foldername/.
	bash RUNME.sh
	probes=( 1a small bencene  )
	for probe in "${probes[@]}"
	do
		mv output-open-BS-SITE-$probe $foldername/.
		mv output-collapsed-BS-SITE-$probe $foldername/.
	done
	unlink jedi.params
done
cp jedi.params.backup jedi.params





