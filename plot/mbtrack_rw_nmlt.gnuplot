set output 'mbtrack_rw_nmlt.eps'
set terminal postscript eps color solid enhanced

# Line style
set style line 1 pt 6 ps 1 lt 1 lw 2 # --- circles points

# Grid style
set style line 12 lc rgb '#808080' lt 0 lw 1
set grid back ls 12

# Borders style
set style line 11 lc rgb '#808080' lt 1
set border 3 back ls 11
set tics nomirror

set title "Growth Time of VERtical amplitude invariant versus Nmlt\n\nSOLEIL: UNIFORM filling RW only (500 mA)"
set xlabel 'Nmlt'
set ylabel 'VERtical growth time (ms)'
plot 'mbtrack_rw_nmlt.dat' using 1:2 with lp ls 1 title 'Center of mass'
