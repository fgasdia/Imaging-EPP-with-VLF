set terminal pngcairo size 700,600 lw 1 fontscale 1
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

# ARG1 = basedir, e.g. letkf_day2/
# ARG2 = filename (including png)
# ARG3 = hlow
# ARG4 = hhi
# ARG5 = blow
# ARG6 = bhi

lm = 0.11
rm = 0.98
bm = 0.13
tm = 0.92
xspace = 0.07
yspace = 0.05
set multiplot layout 3,2 columnsfirst margins lm,rm,bm,tm spacing xspace,yspace

set border 2 front ls 101

set format x '%h'
set format y '%h'

hlim = 5
blim = 0.25

xmax = 8
set xrange [-1:xmax]
set yrange [ARG3:ARG4]  # hlow:hhi

set xtics scale 0
unset xtics

# set label '$\lambda$' at 600,30e3

set label "h' est (km)" at screen 0.11,0.96

numx = 0.03
set label 3 'Grid #' at screen numx,0.96 offset char -2

set label 2 '1' at screen numx,0.81
idx = 0
plot for [i=4:*] ARG1.'hens.csv' u 0:i index idx w p pt 6 ps 0.8 lc rgb '#bb888888' notitle, \
    for [i=1:xmax] ARG1.'hens_gaus.csv' u (column(i+1)+(i-1)):1 index idx w l lc 4 notitle, \

unset label

set label 2 '2' at screen numx,0.54
idx = 1
plot for [i=4:*] ARG1.'hens.csv' u 0:i index idx w p pt 6 ps 0.8 lc rgb '#bb888888' notitle, \
    for [i=1:xmax] ARG1.'hens_gaus.csv' u (column(i+1)+(i-1)):1 index idx w l lc 4 notitle, \

set xlabel 'Iteration'
set xtics 0,1,6 scale 0

set label 2 '3' at screen numx,0.26
idx = 2
plot for [i=4:*] ARG1.'hens.csv' u 0:i index idx w p pt 6 ps 0.8 lc rgb '#bb888888' notitle, \
    for [i=1:xmax] ARG1.'hens_gaus.csv' u (column(i+1)+(i-1)):1 index idx w l lc 4 notitle, \

###
unset xtics
unset xlabel
unset label 2
set yrange [ARG5:ARG6]

# the area-normalized beta gaussians are huge
bscale = 0.03

set label '\beta est (km^{-1})' at screen 0.58,0.96

idx = 0
plot for [i=4:*] ARG1.'bens.csv' u 0:i index idx w p pt 6 ps 0.8 lc rgb '#bb888888' notitle, \
    for [i=1:xmax] ARG1.'bens_gaus.csv' u (column(i+1)*bscale+(i-1)):1 index idx w l lc 1 notitle, \

unset label

idx = 1
plot for [i=4:*] ARG1.'bens.csv' u 0:i index idx w p pt 6 ps 0.8 lc rgb '#bb888888' notitle, \
    for [i=1:xmax] ARG1.'bens_gaus.csv' u (column(i+1)*bscale+(i-1)):1 index idx w l lc 1 notitle, \

set xlabel 'Iteration'
set xtics 0,1,6 scale 0

idx = 2
plot for [i=4:*] ARG1.'bens.csv' u 0:i index idx w p pt 6 ps 0.8 lc rgb '#bb888888' notitle, \
    for [i=1:xmax] ARG1.'bens_gaus.csv' u (column(i+1)*bscale+(i-1)):1 index idx w l lc 1 notitle, \
  
unset multiplot
