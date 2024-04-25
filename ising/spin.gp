#set terminal qt 
set terminal pngcairo
set output 'ising_plot.png'

set title "Average Spin vs Iteration Number"
set xlabel "Iteration Number"
set ylabel "Average Spin"

#set xrange [0:100000]
#set yrange [-1:1]

plot "OUT" using 1:2 with points title "Average Spin"

#pause -1
