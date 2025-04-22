#!/bin/bash
clear

OUTPUT_EGF=false
OUTPUT_NUT=false
OUTPUT_CHEM=false
OUTPUT_OUT=false
OUTPUT_PHENOTYPE=false
OUTPUT_STOCHASTIC=true
OUTPUT_STOCHASTIC_CHEM=true
OUTPUT_ANIMATION=false

echo "Compiling plot source files..."

rm ./build/main_plot.exe

cd src_plot/
g++ PlotFactory.cpp Vector.cpp main.cpp -o ../build/main_plot.exe
cd ../

if ${OUTPUT_EGF}; then
	echo "Plotting EGF files..."

	cd output/egf/

	rm egf*.pdf

	# ./lista.sh
	EGF_FILES=`ls egf*.dat | wc -l`
	for i in $(seq 1 10 ${EGF_FILES})
	do
		INPUT_FILENAME=`ls egf*.dat | tr -s ' ' '\n' | awk NR==${i}`
		OUTPUT_FILENAME=`ls egf*.dat | tr -s ' ' '\n' | awk NR==${i} | cut -d. -f1`
		../../build/main_plot.exe EGF ${INPUT_FILENAME} ${OUTPUT_FILENAME}
		gnuplot ${OUTPUT_FILENAME}.plt
		ps2pdf ${OUTPUT_FILENAME}.eps &> /dev/null
		pdfcrop ${OUTPUT_FILENAME}.pdf ${OUTPUT_FILENAME}.pdf
	done

	rm egf*.plt
	rm egf*.eps

	cd ../../
fi

if ${OUTPUT_NUT}; then
	echo "Plotting NUT files..."

	cd output/nut/

	rm nut*.pdf

	# ./lista.sh
	NUT_FILES=`ls nut*.dat | wc -l`
	for i in $(seq 1 10 ${NUT_FILES})
	do
		INPUT_FILENAME=`ls nut*.dat | tr -s ' ' '\n' | awk NR==${i}`
		OUTPUT_FILENAME=`ls nut*.dat | tr -s ' ' '\n' | awk NR==${i} | cut -d. -f1`
		../../build/main_plot.exe NUT ${INPUT_FILENAME} ${OUTPUT_FILENAME}
		gnuplot ${OUTPUT_FILENAME}.plt
		ps2pdf ${OUTPUT_FILENAME}.eps &> /dev/null
		pdfcrop ${OUTPUT_FILENAME}.pdf ${OUTPUT_FILENAME}.pdf
	done

	rm nut*.plt
	rm nut*.eps

	cd ../../
fi

if ${OUTPUT_CHEM}; then
	echo "Plotting CHEM files..."

	cd output/chem/

	rm chem*.pdf

	# ./lista.sh
	CHEM_FILES=`ls chem*.dat | wc -l`
	for i in $(seq 1 10 ${CHEM_FILES})
	do
		INPUT_FILENAME=`ls chem*.dat | tr -s ' ' '\n' | awk NR==${i}`
		OUTPUT_FILENAME=`ls chem*.dat | tr -s ' ' '\n' | awk NR==${i} | cut -d. -f1`
		../../build/main_plot.exe CHEM ${INPUT_FILENAME} ${OUTPUT_FILENAME}
		gnuplot ${OUTPUT_FILENAME}.plt
		ps2pdf ${OUTPUT_FILENAME}.eps &> /dev/null
		pdfcrop ${OUTPUT_FILENAME}.pdf ${OUTPUT_FILENAME}.pdf
	done

	rm chem*.plt
	rm chem*.eps

	cd ../../
fi

if ${OUTPUT_OUT}; then
	echo "Plotting frame files..."

	cd output/out/

	rm out*.pdf
	rm out*.png

	# ./lista.sh
	OUT_FILES=`ls out*.dat | wc -l`
	for i in $(seq 1 10 ${OUT_FILES})
	do
		INPUT_FILENAME=`ls out*.dat | tr -s ' ' '\n' | awk NR==${i}`
		OUTPUT_FILENAME=`ls out*.dat | tr -s ' ' '\n' | awk NR==${i} | cut -d. -f1`
		../../build/main_plot.exe OUT ${INPUT_FILENAME} ${OUTPUT_FILENAME}
		latex ${OUTPUT_FILENAME}.tex &> /dev/null
		dvips ${OUTPUT_FILENAME}.dvi &> /dev/null
		ps2pdf ${OUTPUT_FILENAME}.ps &> /dev/null
		pdfcrop ${OUTPUT_FILENAME}.pdf ${OUTPUT_FILENAME}.pdf
		pdftoppm ${OUTPUT_FILENAME}.pdf out2-2D-$(( (${i}-1)/10 )) -png -f 1 -singlefile -rx 300 -ry 300
	done

	rm out*.tex
	rm out*.log
	rm out*.aux
	rm out*.dvi
	rm out*.ps

	cd ../../
fi

if ${OUTPUT_PHENOTYPE}; then
	echo "Plotting phenotype file..."

	cd output/out/

	rm -r phen/
	mkdir phen/

	../../build/main_plot.exe PHEN phenotype.dat phenotype
	mv phenotype.dat phen/
	mv phenotype.plt phen/

	cd phen/
	gnuplot phenotype.plt
	pdfcrop phenotype.pdf phenotype.pdf

	rm phenotype.plt

	cd ../../../
fi

if ${OUTPUT_STOCHASTIC}; then
	echo "Plotting stochastic file..."

	cd output/

	rm -r stoc/
	mkdir stoc/

	SEED_FOLDERS=`find seed=* -maxdepth 0 -type d | wc -l`
	for i in $(seq 1 ${SEED_FOLDERS})
	do
		cd seed=${i}/

		../../build/main_plot.exe PHEN phenotype${i}.dat phenotype${i}
		mv phenotype${i}.dat ../stoc/
		rm phenotype${i}.plt

		cd ../
	done

	cd stoc/

	../../build/main_plot.exe STOC stochastic.dat stochastic
	gnuplot stochastic.plt
	pdfcrop stochastic.pdf stochastic.pdf

	# rm stochastic.plt

	cd ../../
fi

if ${OUTPUT_STOCHASTIC_CHEM}; then
	echo "Plotting stochastic CHEM file..."

	cd output/

	rm -r stoc_chem/
	mkdir stoc_chem/

	SEED_FOLDERS=`find seed=* -maxdepth 0 -type d | wc -l`
	for i in $(seq 1 ${SEED_FOLDERS})
	do
		cd seed=${i}/chem/

		../../../build/main_plot.exe CHEM_MEAN chem_mean${i}.dat chem_mean${i}
		mv chem_mean${i}.dat ../../stoc_chem/

		cd ../../
	done

	cd stoc_chem/

	../../build/main_plot.exe STOC_CHEM stochastic_chem.dat stochastic_chem
	gnuplot stochastic_chem.plt
	pdfcrop stochastic_chem.pdf stochastic_chem.pdf

	# rm stochastic_chem.plt

	cd ../../
fi

if ${OUTPUT_ANIMATION}; then
	echo "Generating animation..."

	cd output/out/

	rm -r anim/
	mkdir anim/

	# ./conversao.sh
	OUT_FILES=`ls out*.pdf | wc -l`
	for i in $(seq 1 ${OUT_FILES})
	do
		INPUT_FILENAME=`ls out*.pdf | tr -s ' ' '\n' | awk NR==${i}`
		OUTPUT_FILENAME=`ls out*.pdf | tr -s ' ' '\n' | awk NR==${i} | cut -d. -f1`
		convert -density 500 ${INPUT_FILENAME} ${OUTPUT_FILENAME}.png
		mv ${OUTPUT_FILENAME}.png anim/
	done
	cd anim/
	mencoder 'mf://*.png' -mf type=png:fps=8 -ovc lavc -o agents.avi
	
	rm out*.png

	cd ../../../
fi
