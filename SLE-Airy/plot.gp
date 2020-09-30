# x vs. t
set key off
set xlabel 'x'
set ylabel 'u'
set terminal svg
set out "V544-plot.svg"
set title 'First eigenvector for SLE Airy' font 'Liberation Sans Bold,16'
set style line 1 lc rgb '#0f60ad' lt 1 lw 1 pt 0 ps 1.0
plot 'V544-output.txt' with linespoints ls 1 lt -1
set out
