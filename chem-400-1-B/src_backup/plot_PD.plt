set terminal postscript eps enhanced color font 'Helvetica,10'
set output "output/plot_PD.eps"
set title "PD Curve"
m = "output/PD.dat"
set xlabel "u (nM)"
set ylabel "Drug effect factor"
set autoscale
set xtic auto
set ytic auto
set key off
set border linewidth 2
set style line 1 linecolor rgb "#2196F3" dashtype 3 linewidth 2 pointtype 7 pointsize 0.5
plot m using 1:2 with linespoints linestyle 1
