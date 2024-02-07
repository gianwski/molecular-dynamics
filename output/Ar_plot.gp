set terminal pngcairo enhanced font "arial,15" size 1200,900
set output 'energy_plot_combined.png'

set title "Energies vs Time (Ar)"
set xlabel "Time (fs)"
set ylabel "Energy (hartree)"

plot 'Ar_energy.dat' using 2:3 with lines lw 3 title 'Total Energy', \
     'Ar_energy.dat' using 2:4 with lines lw 3 title 'Kinetic Energy', \
     'Ar_energy.dat' using 2:5 with lines lw 3 title 'Potential Energy'
