set term pdfcairo

set output 'offdiag_ravg.pdf'
set title 'Random access: normalised running average of the singlet yield'
plot 'fort.10' u :1 w lines t 'ravg'

set output 'offdiag_rsum.pdf'
set title 'Random access: running sum of the singlet yield -- Upward curve b/c same trend in nbr of products'
plot 'fort.10' u :2 w lines t 'rsum'

set output 'offdiag_rcount.pdf'
set title 'Random access: running count of element products -- Upward curve is not normal.'
plot 'fort.10' u :3 w lines t 'rcount'

set output 'offdiag_syields.pdf'
set title 'offdiag, parallel access: kernel y value at each step'
plot 'fort.11' u :1 w lines t 'y'

quit
