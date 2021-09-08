set term pdfcairo

set output 'offdiag_ravg.pdf'
set title 'Random access: running average of the singlet yield'
plot 'fort.10' u :1 w lines t 'ravg'

set output 'offdiag_rsum.pdf'
set title 'Random access: running sum of the singlet yield -- Why the upward curve?'
plot 'fort.10' u :2 w lines t 'rsum'

set output 'offdiag_syields.pdf'
set title 'Kernel y value at each step'
plot 'fort.11' u :1 w lines t 'y'

quit
