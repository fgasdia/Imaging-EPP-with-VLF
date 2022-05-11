set terminal pngcairo size 700,600 lw 1 fontscale 1
set output ARG7

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
# load 'amp.pal'
load ARG8

# ARG1 = datadir
#datadir = 'C:\\dev\\SubionosphericVLFInversions\\TypicalVFSA\\daytime_exact_lwpc_progress1\\'
# ARG2 = string label (iteration number or truth)
# ARG3 = column name
# ARG4 = colorbar label
# ARG5 = cb low
# ARG6 = cb high
# ARG7 = filename (including png)
# ARG8 = colormap

mapfile(x) = basemapdir.x."/".x.".gp.csv"

unset xtics
unset ytics
unset border

set lmargin 1
set rmargin 1

set size ratio -1

xmin = -2600000
xmax = 500000
ymin = 500000
ymax = 3300000

set xrange [xmin:xmax]
set yrange [ymin:ymax]
set cbrange [ARG5:ARG6]

set cbtics offset char -1
set cblabel ARG4 offset char -0.5

label_x0 = -400000
label_y0 = 2630000
label_x1 = 250000
label_y1 = 2730000

set multiplot

set label 1 'Iteration: '.ARG2 at label_x0,label_y0 left offset char 0.5,0.5 tc 'black' front
set object 1 rectangle from label_x0,label_y0 to label_x1,label_y1 fillstyle solid noborder fc 'white' front

if ((ARG3 eq "herr") || (ARG3 eq "berr")) {
    plot ARG1.'maps_'.ARG2.'.gp.csv' u 1:2:ARG3 w image notitle, \
        mapfile('ne_50m_graticules_10') w l ls 102 notitle, \
        mapfile('ne_50m_coastline') w l ls 101 notitle, \
        mapfile('ne_50m_admin_1_states_provinces_lakes') w l ls 101 notitle, \
        ARG1.'gcps.gp.csv' w l ls 104 notitle, \
        ARG1.'transmitters.gp.csv' w p pt 9 ps 2 lc 'black' notitle, \
        ARG1.'receivers.gp.csv' w p pt 7 ps 1.2 lc 'black' notitle, \
        ARG1.'ctrlpts_'.ARG2.'.gp.csv' w p pt 6 ps 1 lc 'black' notitle, \
        ARG1.ARG3.'cntr_'.ARG2.'.gp.csv' w l lc 'black' dt 3 lw 1.5 notitle
} else {
    plot ARG1.'maps_'.ARG2.'.gp.csv' u 1:2:ARG3 w image notitle, \
        mapfile('ne_50m_graticules_10') w l ls 102 notitle, \
        mapfile('ne_50m_coastline') w l ls 101 notitle, \
        mapfile('ne_50m_admin_1_states_provinces_lakes') w l ls 101 notitle, \
        ARG1.'gcps.gp.csv' w l ls 104 notitle, \
        ARG1.'transmitters.gp.csv' w p pt 9 ps 2 lc 'black' notitle, \
        ARG1.'receivers.gp.csv' w p pt 7 ps 1.2 lc 'black' notitle, \
        ARG1.'ctrlpts_'.ARG2.'.gp.csv' w p pt 6 ps 1 lc 'black' notitle, \
        ARG1.'cntr_98.gp.csv' w l lc 'white' dt 2 lw 2 notitle
}
    
# plot ARG1.'E.gp.csv' every ::ARG2::ARG2 u (label_x0):(label_y1):(sprintf("E: %.2g", $1)) w labels left offset char 0.8,-2 tc 'black' notitle

unset label 1
unset multiplot
