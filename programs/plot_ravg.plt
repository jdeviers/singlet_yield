set term pdfcairo

set output 'test.pdf'
set title 'Random access: normalised running average of the singlet yield'
plot 'fort.10' u :1 w lines t 'ravg'

set title 'Random access: running sum of the singlet yield'
plot 'fort.10' u :2 w lines t 'rsum'

set title 'Random access: running count of element products'
plot 'fort.10' u :3 w lines t 'rcount'

quit
