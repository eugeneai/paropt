# set terminal postscript portrait enhanced color dashed lw 1 "Helvetica" 14
# set terminal postscript portrait enhanced monochrome dashed lw 1 "Helvetica" 14

set terminal postscript eps enhanced monochrome dashed lw 1 "Helvetica" 24
set output "my-plot.ps"
set termoption lw 4



#set linestyle 1 lt 1 lw 3 pt 5
#set linestyle 2 lt 6 lw 3 pt 6

plot "test_with_plt10.dat" using 1:2 with lines title "X" lt 1 , \
   "test_with_plt10.dat" using 1:3 every 12 with linespoints title "U" lt 6 pt 5 ps 3

#replot
#set terminal x11
#set size 1,1
