#!/bin/bash
clear

CONFIG_FILE="config-2D.cfg"
SEEDS=10
echo "Configuration file = $CONFIG_FILE"
echo "Number of executions = $SEEDS"
echo ""

make clean
make

for SEED in $(seq 21 $SEEDS)
do
	mkdir output/seed=$SEED
	mkdir output/seed=$SEED/chem
	mkdir output/seed=$SEED/egf
	mkdir output/seed=$SEED/nut
done

for SEED in $(seq 21 $SEEDS)
do
	echo ""
	echo "SEED = $SEED"
	echo ""

	./build/main.exe $CONFIG_FILE $SEED

	cd output/
	mv out/*.dat seed=$SEED/
	mv chem/*.dat seed=$SEED/chem/
	mv egf/*.dat seed=$SEED/egf/
	mv nut/*.dat seed=$SEED/nut/
	cd ..
done
