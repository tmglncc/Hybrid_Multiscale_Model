set terminal postscript eps enhanced color font 'Helvetica,10'
set output "output/plot_PK.eps"
set title "PK Curve"
m = "output/PK.dat"
set xlabel "t (seconds)"
set ylabel "u (nM)"
set autoscale
set xtic auto
set ytic auto
set key off
set border linewidth 2
set style line 1 linecolor rgb "#2196F3" dashtype 3 linewidth 2 pointtype 7 pointsize 0.5
# set logscale y
plot m using 1:2 with linespoints linestyle 1
