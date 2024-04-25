#!/bin/bash

directory="${1:-Frames}"   #default to "Frames"

# Loop through numbered frame files
for ((i=1; ; i++)); do
    # Check if frame file exists
    if [ -f "$directory/frame$i.txt" ]; then
        gnuplot -e "filename='$directory/frame$i.txt'" ising.gp
    else
        break
    fi
done

# Create animated GIF
convert -delay 10 -loop 0 "$directory/frame*.png" animated.gif

