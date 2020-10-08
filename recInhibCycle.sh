# cancer cell cycle = 30hr, macrophage recruitment rate = 1e-8
# cycle treatment
set=0
for i in 1 2 4 6 8 10 12 14
do
	for j in 5 10 15 20 25
	do
		./eim baselineTreatments/recruitmentInhibition/cycled $set 100 1e-8 30 0 100 1.0 $j $i
		let "set=set+1"
	done
done

# cancer cell cycle = 20hr, macrophage recruitment rate = 1e-8
# cycle treatment
set=0
for i in 1 2 4 6 8 10 12 14
do
	for j in 5 10 15 20 25
	do
		./eim increasedProlTreatments/recruitmentInhibition/cycled $set 100 1e-8 20 0 100 1.0 $j $i
		let "set=set+1"
	done
done

# cancer cell cycle = 30hr, macrophage recruitment rate = 2e-8
# cycle treatment
set=0
for i in 1 2 4 6 8 10 12 14
do
	for j in 5 10 15 20 25
	do
		./eim increasedRecTreatments/recruitmentInhibition/cycled $set 100 2e-8 30 0 100 1.0 $j $i
		let "set=set+1"
	done
done

# cancer cell cycle = 20hr, macrophage recruitment rate = 2e-8
# cycle treatment
set=0
for i in 1 2 4 6 8 10 12 14
do
	for j in 5 10 15 20 25
	do
		./eim increasedBothTreatments/recruitmentInhibition/cycled $set 100 2e-8 20 0 100 1.0 $j $i
		let "set=set+1"
	done
done