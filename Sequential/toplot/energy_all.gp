reset
set xlabel 'time'
set ylabel 'Energy'
set title 'Different Energies VS. Time'
filename(n) = sprintf("data.out_%dpart",n)

plot for [i=4:12:2] filename(i) using 1:2 with lines title 'Potential Energy'
pause mouse any
plot for [i=4:12:2] filename(i) using 1:3 with lines title 'Kinetic Energy'
pause mouse any
plot for [i=4:12:2] filename(i) using 1:4 with lines title 'Mechanical Energy'
pause mouse any
reset
