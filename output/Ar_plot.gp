set terminal pngcairo enhanced font "arial,10" size 800,600
set output 'energy_plot_combined.png'

set title "Energies vs Time (Ar)"
set xlabel "Time (fs)"
set ylabel "Energy (hartree)"

plot 'energy' using 2:3 with lines title 'Total Energy', \
     'energy' using 2:4 with lines title 'Kinetic Energy', \
     'energy' using 2:5 with lines title 'Potential Energy'
