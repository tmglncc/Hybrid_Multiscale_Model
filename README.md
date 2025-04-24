# Hybrid Multiscale Model

Cancer is a group of diseases characterized by complex phenomena across multiple temporal and spatial scales. Comprehending its growth dynamics is a challenge and may improve the understanding of underlying mechanisms, and suggests more effective therapy protocols. In this context, mathematical and computational modeling may provide insights into tumorigenesis, cancer growth, and response to drug treatments. In this work, we develop a hybrid discrete-continuum model describing the avascular phase of cancer growth and incorporate chemotherapeutic drugs acting in different phases of the cell cycle. The growth phenomenon occurs at three scales: (i) at the tissue scale, partial differential equations model oxygen, drug, and growth factor dispersion; (ii) at the cellular scale, an agent-based model describes transitions among phenotypic states of each tumor cell and mechanical interactions among cells and the microenvironment; (iii) at the molecular level, ordinary differential equations represent signaling pathways that regulate cellular metabolism, cell cycle, and cell proliferation. Computational experiments demonstrate that the proposed modeling framework can be instrumental in the development of innovative new treatments for cancer patients.

![Schematic representation of the hybrid multiscale model developed](https://drive.google.com/uc?export=view&id=1eaA7yYbQnCQ8Qp6K1iNr5KlRQgk5UxGo)

## Requirements

Our experiments have been performed using **C++ 7.5.0**, **Shell Script 5.0.16**, and **Gnuplot 5.2**.

## Running experiments

1. Clone this repository directly from terminal:
	 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`$ git clone https://github.com/tmglncc/Hybrid_Multiscale_Model.git`
	
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**OR**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Download the _.zip_ file and decompress it.

2. Enter the directory of the experiment you would like to run:
   - _chem-400-1-B_: experiment with application of cisplatin at _t = 400_ h;
   - _chem-400-1-C_: experiment with application of taxotere at _t = 400_ h; or
   - _no-chem_: control experiment.

3. If necessary, make changes in the hybrid multiscale model parameters by modifying the _config-2D.cfg_ file. Set the number of replicates to be executed in the _run_all.sh_ file.

4. Clean the experiment's object files and executable file:
	 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`$ make clean`

5. Compile the experiment's source codes:
	 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`$ make`

6. Run one replicate for the experiment:
	
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`$ make run`
	
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**OR**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Run multiple replicates for the experiment:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`$ ./run_all.sh`

7. If multiple replicates were executed:
    1. Set the types of plots to be generated in the _plot_all.sh_ file. If necessary, make changes in the plot parameters by modifying the _src_plot/PlotFactory.hpp_ file.
    2. Run the script to generate the plots:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`$ ./plot_all.sh`

8. Enter the experiment's _output_ folder and check out the results.

## Cite as

NAOZUKA, G. T.; ROCHA, H. L.; ALMEIDA, R. C. Hybrid multiscale modeling of tumor growth with chemotherapeutic drug dispersion. In: _XLI Ibero-Latin-American Congress on Computational Methods in Engineering (CILAMCE)_. Foz do Iguaçu, PR: ABMEC, 2020. p. 1–7. ISSN 2675-6269. Available in: [https://publicacoes.softaliza.com.br/cilamce2020/article/view/6425](https://publicacoes.softaliza.com.br/cilamce2020/article/view/6425)

Naozuka, G.T.; Rocha, H.L.; Almeida, R.C. Hybrid Multiscale Model, 2020. Version 1.0. Available online: [https://github.com/tmglncc/Hybrid_Multiscale_Model](https://github.com/tmglncc/Hybrid_Multiscale_Model) (accessed on 24 April 2025), doi: [10.5281/zenodo.15276349](https://doi.org/10.5281/zenodo.15276349)
