# y vs. x
set key off
set xlabel 'x'
set ylabel 'y'
set terminal svg
set out "RandomExample2.svg"
set title 'RandomExample2; k=1' font 'Liberation Sans Bold,16'
set style line 1 lc rgb '#0060ad' lt 2 lw 2 pt 0 ps 1.0
plot 'Sol-y.txt' with linespoints ls 1
set out

# dy/dx vs x
set key off
set xlabel 'x'
set ylabel 'dy/dx'
set terminal svg
set out "RandomExample2-dy.svg"
set title 'RandomExample2 dy/dx' font 'Liberation Sans Bold,16'
set style line 1 lc rgb '#0060ad' lt 2 lw 2 pt 0 ps 1.0
plot 'Sol-dy.txt' with linespoints ls 1
set out

# dy/dx vs y
set key off
set xlabel 'y'
set ylabel 'dy/dx'
set terminal svg
set out "RandomExample2-phase.svg"
set title 'RandomExample2 phase' font 'Liberation Sans Bold,16'
set style line 1 lc rgb '#0060ad' lt 2 lw 2 pt 0 ps 1.0
plot 'Sol-phase.txt' with linespoints ls 1
set out
