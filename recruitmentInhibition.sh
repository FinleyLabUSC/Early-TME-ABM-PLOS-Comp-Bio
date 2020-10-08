# cancer cell cycle = 30hr, macrophage recruitment rate = 1e-8
# continuous treatment
set=0
for i in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
do
	./eim baselineTreatments/recruitmentInhibition/continuous $set 100 1e-8 30 0 100 $i 0 0
	let "set=set+1"
done

# cancer cell cycle = 20hr, macrophage recruitment rate = 1e-8
# continuous treatment
set=0
for i in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
do
	./eim increasedProlTreatments/recruitmentInhibition/continuous $set 100 1e-8 20 0 100 $i 0 0
	let "set=set+1"
done

# cancer cell cycle = 30hr, macrophage recruitment rate = 2e-8
# continuous treatment
set=0
for i in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
do
	./eim increasedRecTreatments/recruitmentInhibition/continuous $set 100 2e-8 30 0 100 $i 0 0
	let "set=set+1"
done

# cancer cell cycle = 20hr, macrophage recruitment rate = 2e-8
# continuous treatment
set=0
for i in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
do
	./eim increasedBothTreatments/recruitmentInhibition/continuous $set 100 2e-8 20 0 100 $i 0 0
	let "set=set+1"
done