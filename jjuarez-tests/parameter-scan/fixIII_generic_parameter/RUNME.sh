#Probes
probes=( 1a small bencene  )

for probe in "${probes[@]}"
do
	#Clean up
	rm -r output-open-BS-SITE-$probe
	rm -r output-collapsed-BS-SITE-$probe
	mkdir output-open-BS-SITE-$probe
	mkdir output-collapsed-BS-SITE-$probe

	#========= WORKING ON THE OPEN  BINDING SITE STRUCTURE ========#
	cd output-open-BS-SITE-$probe
	cp ../input/* .
	cp ../probes/* .
	rm simulation2_cluster14800_aligned.pdb

	# Create the grid the polar and apolar atoms for initial structure
	../jedi-setup.py -i initial.pdb -l ligand.pdb -c 0.6 -s 0.15 -a apolar_initial.pdb -p polar_initial.pdb -g grid_initial.pdb -ln LIG
	
	#select the probe file in the PLUMED input
	sed "s/XXXX/probe_$probe.pdb/g" driver_initial-SITE.dat > temp
	mv temp driver_initial-SITE.dat  
	# Caluclate the JEDI score for the initial.pdb conformation using files obtained from the inital structure
	plumed driver --mf_pdb initial.pdb --plumed=driver_initial-SITE.dat

	python ../grid_analysis.py -i grid-step-0.xyz -a acti-step-0.txt -o activity-grid-initial-SITE.pdb

	cd ..
	#========= WORKING ON THE COLLAPSED BINDING SITE STRUCTURE ========#
	cd output-collapsed-BS-SITE-$probe

	cp ../input/* .
	cp ../probes/* .
	rm initial.pdb

	# Create the grid the polar and apolar atoms for 14800 structure
	../jedi-setup.py -i simulation2_cluster14800_aligned.pdb -l ligand.pdb -c 0.6 -s 0.15 -a apolar_14800.pdb -p polar_14800.pdb -g grid_14800.pdb -ln LIG

	#select the probe file in the PLUMED input
        sed "s/XXXX/probe_$probe.pdb/g" driver_14800-SITE.dat > temp
        mv temp driver_14800-SITE.dat

	# Calculate the JEDI score for the initial.pdb conformation using files obtained from the 14800 structure
	plumed driver --mf_pdb simulation2_cluster14800_aligned.pdb --plumed=driver_14800-SITE.dat

	python ../grid_analysis.py -i grid-step-0.xyz -a acti-step-0.txt -o activity-grid-14800.pdb

	cd ..
done


