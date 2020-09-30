# x vs. t
set key off
set xlabel 't'
set ylabel 'x'
set terminal svg
set out "randomexample1-standard.svg"
set title 'Random Example 1' font 'Liberation Sans Bold,16'
set style line 1 lc rgb '#0060ad' lt 2 lw 2 pt 0 ps 1.0
plot 'randomexample1.txt' with linespoints ls 1
set out

# dx/dt vs t
set key off
set xlabel 't'
set ylabel 'dx/dt'
set terminal svg
set out "randomexample1-derivative.svg"
set title 'Random Example 1' font 'Liberation Sans Bold,16'
set style line 1 lc rgb '#0060ad' lt 2 lw 2 pt 0 ps 1.0
plot 'randomexample1-dx.txt' with linespoints ls 1
set out

# dx/dt vs x
set key off
set xlabel 'x'
set ylabel 'dx/dt'
set terminal svg
set out "randomexample1-phase.svg"
set title 'Random Example 1' font 'Liberation Sans Bold,16'
set style line 1 lc rgb '#0060ad' lt 2 lw 2 pt 0 ps 1.0
plot 'randomexample1-phase.txt' with linespoints ls 1
set out
