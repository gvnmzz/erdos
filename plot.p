set term pdf enhanced color font "Ubuntu,10"

set style line 1 linecolor rgb '#e6e6e6' lw 2 ps 0.5
set style line 2 linecolor rgb '#3abbc9' lw 2 ps 0.5
set style line 3 linecolor rgb '#9bca3e' lw 2 ps 0.5
set style line 4 linecolor rgb '#feeb51' lw 2 ps 0.5
set style line 5 linecolor rgb '#ffb92a' lw 2 ps 0.5
set style line 6 linecolor rgb '#ed5314' lw 2 ps 0.5

set key top right spacing 1
set size square
set xlabel 'p'
set ylabel 'Fraction of Singular Graphs'
set output 'fraction.pdf'

plot 'n10.res' u 1:2 w lp ls 6 t 'n=10',\
'n100.res' u 1:2 w lp ls 5 t 'n=100',\
'n300.res' u 1:2 w lp ls 3 t 'n=300',\
'n1000.res' u 1:2 w lp ls 2 t 'n=100'

set output

set xlabel 'pn^{5/8}/log{(n)}'
set ylabel 'Fraction of Singular Graphs'
set output 'fraction_scaled.pdf'

plot 'n10.res' u ($1*10**(5./8.)/log(10)):2 w lp ls 6 t 'n=10',\
'n100.res' u ($1*100**(5./8.)/log(100)):2 w lp ls 5 t 'n=100',\
'n300.res' u ($1*300**(5./8.)/log(300)):2 w lp ls 3 t 'n=300',\
'n1000.res' u ($1*1000**(5./8.)/log(1000)):2 w lp ls 2 t 'n=100'

set output
