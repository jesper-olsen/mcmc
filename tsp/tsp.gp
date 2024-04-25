set term qt
set xlabel 'X'
set ylabel 'Y'

set title 'Traveling Salesman Problem'
plot 'cities.dat' u 2:1:(sprintf("%d", $0+1)) with labels point pt 7 offset char 1,1 notitle, \
     'cities.dat' u 2:1 with points pointtype 7 pointsize 1 notitle, \
     'route.dat' u 2:1:(sprintf("%d", $0)) with lines lw 2 title 'Route'

pause -1
