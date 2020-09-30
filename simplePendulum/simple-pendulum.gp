# theta vs. t
set key off
set xlabel 't'
set ylabel 'θ'
set terminal svg
set out "simple-pendulum-standard.svg"
set title 'Simple Pendulum' font 'Liberation Sans Bold,16'
set style line 1 lc rgb '#0060ad' lt 2 lw 2 pt 0 ps 1.0
plot 'simple-pendulum.txt' with linespoints ls 1
set out

# dtheta/dt vs t
set key off
set xlabel 't'
set ylabel 'dθ/dt'
set terminal svg
set out "simple-pendulum-derivative.svg"
set title 'Simple Pendulum' font 'Liberation Sans Bold,16'
set style line 1 lc rgb '#0060ad' lt 2 lw 2 pt 0 ps 1.0
plot 'simple-pendulum-dx.txt' with linespoints ls 1
set out

# dtheta/dt vs theta
set key off
set xlabel 'theta'
set ylabel 'd theta/dt'
set terminal svg
set out "simple-pendulum-phase.svg"
set title 'Simple Pendulum' font 'Liberation Sans Bold,16'
set style line 1 lc rgb '#0060ad' lt 2 lw 2 pt 0 ps 1.0
plot 'simple-pendulum-phase.txt' with linespoints ls 1
set out
