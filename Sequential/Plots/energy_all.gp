reset
set xlabel 'time'
set ylabel 'Energy'
set title 'Different Energies VS. Time'
filename(n) = sprintf("data.out_%dpart",n)

plot for [i=4:12:2] filename(i) using 1:2 with lines title 'Potential Energy'
set term 'png'
set output 'Potential_energy_allpart.png'
replot
unset term
pause mouse any
plot for [i=4:12:2] filename(i) using 1:3 with lines title 'Kinetic Energy'
set term 'png'
set output 'Kinetic_energy_allpart.png'
replot
unset term
pause mouse any
plot for [i=4:12:2] filename(i) using 1:4 with lines title 'Mechanical Energy'
set term 'png'
set output 'Mechanical_energy_allpart.png'
replot
unset term
pause mouse any
reset
