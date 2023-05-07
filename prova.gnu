set term png
set output "prova.png"

file = 'den.dat'

plot file index 0 using 1:2 with lines , \
    file index 1 using 1:2 with lines, \
    file index 2 using 1:2 with lines , \
    file index 3 using 1:2 with lines , \