# cancer cell cycle = 30hr, macrophage recruitment rate = 1e-8
# cycle treatment
set=0
for i in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
do
	for j in 5 10 15 20 25
	do
		./eim baselineTreatments/macrophageDepletion/cycled $set 100 1e-8 30 1 100 $i $j 0.02083333
		let "set=set+1"
	done
done

# cancer cell cycle = 20hr, macrophage recruitment rate = 1e-8
# cycle treatment
set=0
for i in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
do
	for j in 5 10 15 20 25
	do
		./eim increasedProlTreatments/macrophageDepletion/cycled $set 100 1e-8 20 1 100 $i $j 0.02083333
		let "set=set+1"
	done
done

# cancer cell cycle = 30hr, macrophage recruitment rate = 2e-8
# cycle treatment
set=0
for i in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
do
	for j in 5 10 15 20 25
	do
		./eim increasedRecTreatments/macrophageDepletion/cycled $set 100 2e-8 30 1 100 $i $j 0.02083333
		let "set=set+1"
	done
done

# cancer cell cycle = 20hr, macrophage recruitment rate = 2e-8
# cycle treatment
set=0
for i in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
do
	for j in 5 10 15 20 25
	do
		./eim increasedBothTreatments/macrophageDepletion/cycled $set 100 2e-8 20 1 100 $i $j 0.02083333
		let "set=set+1"
	done
done