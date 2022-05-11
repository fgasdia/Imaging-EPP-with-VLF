set terminal pngcairo size 600,300 lw 1 fontscale 1
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

# ARG1 = base dir
# ARG2 = output file
# ARG3 = amp_resid.csv column for last iteration (usually "6" for LETKF)

lm = 0.1
rm = 0.99
bm = 0.04
tm = 0.9
xspace = 0.12
set multiplot

set lmargin at screen lm
set rmargin at screen 0.5
set tmargin at screen tm
set bmargin at screen bm

set border 2 front ls 101

set format x '%h'
set format y '%h'

alim = 20
plim = 180

xmax = ARG3-0.5
set xrange [-0.5:xmax]

posfilt(x) = x <= 0 ? NaN : x
negfilt(x) = x >= 0 ? NaN : x


base = 10.0
linthresh = 0.5
linscale = 0.5

linscale_adj = linscale / (1.0 - base**-1)
log_base = log(base)

# https://matplotlib.org/stable/_modules/matplotlib/scale.html#SymmetricalLogScale
# assumes x is positive
loglin(x) = x <= linthresh ? x*linscale_adj : linthresh*(linscale_adj + log(x/linthresh)/log_base)

set yrange [-loglin(alim):loglin(alim)]

set xtics scale 0
unset xtics
unset ytics  # we do it below so it isn't double plotted


s = loglin(0.1)

# positive
plot (-s) w filledcurves y=s fc rgb '#404040' fs transparent solid 0.5 notitle, \
    (-2*s) w filledcurves y=2*s fc rgb '#a0a0a0' fs transparent solid 0.5 notitle, \
    for [i = 1:0+ARG3] ARG1.'amp_resid.csv' u (i-1):((loglin(posfilt(column(i))))) w p pt 6 lc 1 notitle


# Labels here so we don't need to remove them (from being double printed)
set ytics ("" -loglin(16) 1, "" -loglin(18) 1, "-20" -loglin(20), \
    "-10" -loglin(10), "" -loglin(12) 1, "" -loglin(14) 1, \
    "-5" -loglin(5), "" -loglin(6) 1, "" -loglin(7) 1, "" -loglin(8) 1, "" -loglin(9) 1, \
    "-1" -loglin(1.0), "" -loglin(2) 1, "" -loglin(3) 1, "" -loglin(4) 1, \
    "-0.5" -loglin(0.5), \
    "0" loglin(0), \
    "0.5" loglin(0.5), \
    "1" loglin(1.0), "" loglin(2) 1, "" loglin(3) 1, "" loglin(4) 1, \
    "5" loglin(5), "" loglin(6) 1, "" loglin(7) 1, "" loglin(8) 1, "" loglin(9) 1, \
    "10" loglin(10), "" loglin(12) 1, "" loglin(14) 1, \
    "" loglin(16) 1, "" loglin(18) 1, "20" loglin(20)) nomirror
# 1 means minor tick

set arrow from -0.5,0 to xmax,0 as 101 dt 2 front
set label 1 "H(x) - y_o" at 0.6,loglin(alim) center offset char 0,1

set ylabel 'Amplitude (dB)' offset char 1

set arrow 20 from -0.5,loglin(linthresh) to xmax,loglin(linthresh) as 101 dt 3 back
set arrow 21 from -0.5,-loglin(linthresh) to xmax,-loglin(linthresh) as 101 dt 3 back

# negative
plot for [i = 1:0+ARG3] ARG1.'amp_resid.csv' u (i-1):(-loglin(abs(negfilt(column(i))))) w p pt 6 lc 1 notitle

########
# Phase

set lmargin at screen 0.6
set rmargin at screen rm

unset ytics  # we do it below
unset ylabel
unset for [i=10:11] arrow i
unset label

base = 10.0
linthresh = 5.0
linscale = 0.5

linscale_adj = linscale / (1.0 - base**-1)
log_base = log(base)

# assumes x is positive
loglin(x) = x <= linthresh ? x*linscale_adj : linthresh*(linscale_adj + log(x/linthresh)/log_base)

set yrange [-loglin(plim):loglin(plim)]

s = loglin(1.0)
plot (-s) w filledcurves y=s fc rgb '#404040' fs transparent solid 0.5 notitle, \
    (-2*s) w filledcurves y=2*s fc rgb '#a0a0a0' fs transparent solid 0.5 notitle, \
    for [i = 1:0+ARG3] ARG1.'phase_resid.csv' u (i-1):(loglin(posfilt(column(i)))) w p pt 6 lc 1 notitle
    # for [i = 1:2] 'cobyla_day2/phase_resid.csv' u (i-0.8):(loglin(posfilt(column(cobyla[i])))) w p pt 6 lc 3 notitle


# Labels here so we don't need to remove them (from being double printed)
set ytics ("-180" -loglin(180), \
    "-90" -loglin(90), "" -loglin(120) 1, "" -loglin(150) 1, \
    "-45" -loglin(45), "" -loglin(60) 1, "" -loglin(75) 1, \
    "-10" -loglin(10), "" -loglin(20) 1, "" -loglin(30) 1, \
    "-5" -loglin(5), \
    "0" loglin(0), \
    "5" loglin(5), \
    "10" loglin(10), "" loglin(20) 1, "" loglin(30) 1, \
    "45" loglin(45), "" loglin(60) 1, "" loglin(75) 1, \
    "90" loglin(90), "" loglin(120) 1, "" loglin(150) 1,\
    "180" loglin(180)) nomirror

set arrow from -0.5,0 to xmax,0 as 101 dt 2 front
set label 1 "H(x) - y_o" at 0.6,loglin(plim) center offset char 0,1

set ylabel 'Phase (deg)' offset char 1.5

set arrow 20 from -0.5,loglin(linthresh) to xmax,loglin(linthresh) as 101 dt 3 back
set arrow 21 from -0.5,-loglin(linthresh) to xmax,-loglin(linthresh) as 101 dt 3 back

plot for [i = 1:ARG3] ARG1.'phase_resid.csv' u (i-1):(-loglin(abs(negfilt(column(i))))) w p pt 6 lc 1 notitle

unset multiplot