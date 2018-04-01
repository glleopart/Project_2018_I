reset
set xlabel 'time'
set ylabel 'Temperature'
set title 'Temperature VS. Time'
unset key

filename(n) = sprintf("traj.xyz_%dpart",n)

plot for [i=4:12:2] filename(i) using 1:6 with lines title 'Temperature'
pause mouse any

reset
