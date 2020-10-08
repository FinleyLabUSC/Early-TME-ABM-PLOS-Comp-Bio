# cancer cell cycle = 30hr, macrophage recruitment rate = 1e-8
# continuous treatment
set=0
for i in 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.01
do
	./eim baselineTreatments/macrophageDepletion/continuous $set 100 1e-8 30 1 100 $i 0 0
	let "set=set+1"
done

# cancer cell cycle = 20hr, macrophage recruitment rate = 1e-8
# continuous treatment
set=0
for i in 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.01
do
	./eim increasedProlTreatments/macrophageDepletion/continuous $set 100 1e-8 20 1 100 $i 0 0
	let "set=set+1"
done

# cancer cell cycle = 30hr, macrophage recruitment rate = 2e-8
# continuous treatment
set=0
for i in 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.01
do
	./eim increasedRecTreatments/macrophageDepletion/continuous $set 100 2e-8 30 1 100 $i 0 0
	let "set=set+1"
done

# cancer cell cycle = 20hr, macrophage recruitment rate = 2e-8
# continuous treatment
set=0
for i in 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.01
do
	./eim increasedBothTreatments/macrophageDepletion/continuous $set 100 2e-8 20 1 100 $i 0 0
	let "set=set+1"
done