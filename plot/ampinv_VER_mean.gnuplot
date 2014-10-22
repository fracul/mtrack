set output 'ampinv_VER_mean.eps'
set terminal postscript eps color solid enhanced size 4.0,5.5

# Same as:
#   set multiplot layout 2,1 title "_TITLE_"
# Where _TITLE_ is obtained from a shell command line

load '< echo "set multiplot layout 2,1 title \"VERtical amplitude invariant\\n\\n`grep -o -E "[ ]*title[ ]*[=:][ ]*[^[:cntrl:]]*" *.conf`\\nData file created/modified: `date -r ampinv_VER_mean.dat +%F`\" "'

set title 'Amplitude invariant versus number of turns'
set xlabel 'Turn'
set ylabel 'ampinv'
plot 'ampinv_VER_mean.dat' using 1:2 title 'Bunch center of mass', 'ampinv_VER_mean.dat' using 1:3 title 'Bunch slices'

set title 'Logscale Amplitude invariant versus number of turns'
set logscale y
set xlabel 'Turn'
set ylabel 'log(ampinv)'
set key right bottom
plot 'ampinv_VER_mean.dat' using 1:2 title 'Bunch center of mass', 'ampinv_VER_mean.dat' using 1:3 title 'Bunch slices'

unset multiplot
