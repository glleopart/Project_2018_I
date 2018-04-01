reset
set xlabel 'time'
set ylabel 'Momentum'
unset key
set multiplot layout 3,2 
set tmargin 2
#
set  title 'Momentum Behaviour N=64 particles'
plot 'data.out_4part' using 1:5 with lines
#
set  title 'Momentum Behaviour N=216 particles'
plot 'data.out_6part' using 1:5 with lines
# 
set  title 'Momentum Behaviour N=512 particles'
plot 'data.out_8part' using 1:5 with lines
# 
set  title 'Momentum Behaviour N=1000 particles'
plot 'data.out_10part' using 1:5 with lines
# 
set  title 'Momentum Behaviour N=1728 particles'
plot 'data.out_12part' using 1:5 with lines
# 
unset multiplot
pause mouse any
reset
