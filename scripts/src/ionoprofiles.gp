set terminal pngcairo size 900,500 fontscale 1
set output ARG2

if (GPVAL_SYSNAME eq "Linux") {
    set loadpath '~/gnuplot'
    basemapdir = '/home/forrest/gis/basemaps/'
    load 'src/common.cfg'
} else {
    set loadpath 'C:\\Users\\forrest\\gnuplot\\'
    # set loadpath 'C:\\Users\\lair\\gnuplot\\'
    basemapdir = 'C:\\gis\\basemaps\\'
    load 'src\\common.cfg'
}

load 'mydistinctcolors.pal'

set palette defined ( 0 '#000075', 1 '#ffffff') negative

lm = 0.08
rm = 0.98
bm = 0.12
tm = 0.96
space = 0.02
set multiplot layout 1,3 margins lm,rm,bm,tm spacing space

set xrange [1e3:3*1e11]
set yrange [40:112]

set logscale x
# set xtics 100
unset mxtics

set xlabel 'N_e (m^{-3})'
set format x '10^{%T}'

set format y '%h'
set ylabel 'Altitude (km)'

set grid back ls 102
# set grid back lc rgb '#eeeeee' lw 0.5

# Plot 1 (upper left)
# set label 3 'Day' at 1e4,105 front center boxed bs 2

set style line 11 lt 1 lw 2 
set style line 12 lt 2 lw 2
set style line 13 lt 2 lw 2 dt 2

set title '1' offset 0,-1
plot ARG1.'ionoprofiles.csv' u 'true1':'alt' w l ls 11 notitle, \
	'' u (column('alt') <= column('Wr1') ? column('wait1') : NaN):'alt' w l ls 12 notitle, \
	'' u (column('alt') > column('Wr1') ? column('wait1') : NaN):'alt' w l ls 13 notitle

set ytics scale 0
set ytics format ''
unset ylabel
set border 1

set title '2'
plot ARG1.'ionoprofiles.csv' u 'true2':'alt' w l ls 11 notitle, \
	'' u (column('alt') <= column('Wr2') ? column('wait2') : NaN):'alt' w l ls 12  notitle, \
	'' u (column('alt') > column('Wr2') ? column('wait2') : NaN):'alt' w l ls 13 notitle

set title '3'
plot ARG1.'ionoprofiles.csv' u 'true3':'alt' w l ls 11 notitle, \
	'' u (column('alt') <= column('Wr3') ? column('wait3') : NaN):'alt' w l ls 12  notitle, \
	'' u (column('alt') > column('Wr3') ? column('wait3') : NaN):'alt' w l ls 13 notitle

unset multiplot