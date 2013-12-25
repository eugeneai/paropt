
plot "test_with_plt0.dat" using 1:2 with lines title "X1", \
     "test_with_plt0.dat" using 1:3 with lines title "U1", \
     "test_with_plt1.dat" using 1:2 with lines title "X2", \
     "test_with_plt1.dat" using 1:3 with lines title "U2"
set size 1.0, 0.6
set terminal postscript portrait enhanced mono dashed lw 1 "Helvetica" 8
set output "my-plot.ps"
replot
set terminal x11
set size 1,1
