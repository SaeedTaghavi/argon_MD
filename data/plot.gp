set terminal png
set output 'plot.png'
set grid lw 3
set key off
set xrange [0:20]
set yrange [0:3]
set xlabel '$r$'
set ylabel '$g(r)$'
f(x)=a*x+b
plot 'plot' using 1:2 w lines lw 2
