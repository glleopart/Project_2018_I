reset
set term pngcairo dashed enhanced
set output 'temperature_allpart.png'
set xlabel 'time'
set ylabel 'Temperature'
set title 'Temperature VS. Time'
unset key

set multiplot layout 3,2 title 
set tmargin 2
#
set  title 'N=64 particles'
set yrange [-30:30]
set ytics -30,15,30
plot 'data.out_4part' using 1:6 with lines
#
set  title 'N=216 particles'
set yrange [-60:60]
set ytics -60,30,60
plot 'data.out_6part' using 1:6 with lines
# 
set  title 'N=512 particles'
set yrange [-100:100]
set ytics -100,50,100
plot 'data.out_8part' using 1:6 with lines
# 
set  title ' N=1000 particles'

set yrange [-100:150]
set ytics -100,50,150
plot 'data.out_10part' using 1:6 with lines
# 
set  title ' N=1728 particles'
set yrange [-150:150]
set ytics -150,100,150
plot 'data.out_12part' using 1:6 with lines
#
unset multiplot
pause mouse any
unset term
reset
