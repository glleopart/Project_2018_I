reset
set xlabel 'time'
set ylabel 'Temperature'
set title 'Temperature VS. Time'
unset key
set multiplot layout 3,2 title 'Temperature vs time'
set tmargin 2
#
set  title 'Temperature Behaviour N=64 particles'
plot 'data.out_4part' using 1:6 with lines
#
set  title 'Temperature Behaviour N=216 particles'
plot 'data.out_6part' using 1:6 with lines
# 
set  title 'Temperature Behaviour N=512 particles'
plot 'data.out_8part' using 1:6 with lines
# 
set  title 'Temperature Behaviour N=1000 particles'
plot 'data.out_10part' using 1:6 with lines
# 
set  title 'Temperature Behaviour N=1728 particles'
plot 'data.out_12part' using 1:6 with lines
# 
unset multiplot
pause mouse any
reset
