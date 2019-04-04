for i in `cat params2scan` 
do 
	echo "===================Scanning JEDI parameter: $i==================================="
	cp -r fixIII_generic_parameter fixIII_scan_$i 
	cd fixIII_scan_$i 
	bash RUNSCAN.sh $i
	cd ..

done
