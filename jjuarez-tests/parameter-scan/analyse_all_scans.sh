for i in `cat params2scan` 
do 
	echo "=================== Analysis of JEDI parameter: $i==================================="
	python analyse_scan.py	-f fixIII_scan_$i -p $i -th 1.0
done
