#set terminal qt 
set terminal pngcairo
set output 'ising_plot.png'

# Set the title and labels for the axes
set title "Average Spin vs Iteration Number"
set xlabel "Iteration Number"
set ylabel "Average Spin"

# Set the range for the axes
#set xrange [0:100000]
#set yrange [-1:1]

# Plot the data using lines
plot "OUT" using 1:2 with points title "Average Spin"

#pause -1
