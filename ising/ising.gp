set terminal pngcairo enhanced font 'arial,10' size 800, 650

set title 'Ising Model Spin Configuration'
set xlabel 'Column'
set ylabel 'Row'

set palette defined (-1 "blue", 0 "white", 1 "red")

change_extension(filename, ext) = strcol(filename) . "." . ext

set xrange [0:749]
set yrange [0:599]

set output filename . ".png"
plot filename matrix with image 
