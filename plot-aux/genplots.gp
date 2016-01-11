reset
set term post eps enhanced color dashed "Helvetica" 24
set output 'pwg-NLO.eps-_sigma_total_index___0.eps'
iindex = 0
mytitle = "sigma_total"
infil1 = "pwg-NLO.top"

#set log x
set log y

set datafile fortran
set style data xyerrorlines

set xtics format ""

set xrange [] writeback
set yrange [] 
set ylabel '{d {/Symbol s}}/{d {p_T}} [fb]'
set format y "10^{%L}"

plot infil1 index  iindex using ($1+$2)/2:3:1:2:($3-$4):($3+$4) title mytitle

reset
set term post eps enhanced color dashed "Helvetica" 24
set output 'pwg-NLO.eps-_inv_mass_neutralino_pair_index___1.eps'
iindex = 1
mytitle = "inv_mass_neutralino_pair"
infil1 = "pwg-NLO.top"

#set log x
set log y

set datafile fortran
set style data xyerrorlines

set xtics format ""

set xrange [] writeback
set yrange [] 
set ylabel '{d {/Symbol s}}/{d {p_T}} [fb]'
set format y "10^{%L}"

plot infil1 index  iindex using ($1+$2)/2:3:1:2:($3-$4):($3+$4) title mytitle

reset
set term post eps enhanced color dashed "Helvetica" 24
set output 'pwg-NLO.eps-_pt_neu1_index___2.eps'
iindex = 2
mytitle = "pt_neu1"
infil1 = "pwg-NLO.top"

#set log x
set log y

set datafile fortran
set style data xyerrorlines

set xtics format ""

set xrange [] writeback
set yrange [] 
set ylabel '{d {/Symbol s}}/{d {p_T}} [fb]'
set format y "10^{%L}"

plot infil1 index  iindex using ($1+$2)/2:3:1:2:($3-$4):($3+$4) title mytitle

reset
set term post eps enhanced color dashed "Helvetica" 24
set output 'pwg-NLO.eps-_eta_neu1_index___3.eps'
iindex = 3
mytitle = "eta_neu1"
infil1 = "pwg-NLO.top"

#set log x
set log y

set datafile fortran
set style data xyerrorlines

set xtics format ""

set xrange [] writeback
set yrange [] 
set ylabel '{d {/Symbol s}}/{d {p_T}} [fb]'
set format y "10^{%L}"

plot infil1 index  iindex using ($1+$2)/2:3:1:2:($3-$4):($3+$4) title mytitle

reset
set term post eps enhanced color dashed "Helvetica" 24
set output 'pwg-NLO.eps-_number_of_visible_jets_index___4.eps'
iindex = 4
mytitle = "number_of_visible_jets"
infil1 = "pwg-NLO.top"

#set log x
set log y

set datafile fortran
set style data xyerrorlines

set xtics format ""

set xrange [] writeback
set yrange [] 
set ylabel '{d {/Symbol s}}/{d {p_T}} [fb]'
set format y "10^{%L}"

plot infil1 index  iindex using ($1+$2)/2:3:1:2:($3-$4):($3+$4) title mytitle

reset
set term post eps enhanced color dashed "Helvetica" 24
set output 'pwg-NLO.eps-_pt_hardest_jet_index___5.eps'
iindex = 5
mytitle = "pt_hardest_jet"
infil1 = "pwg-NLO.top"

#set log x
set log y

set datafile fortran
set style data xyerrorlines

set xtics format ""

set xrange [] writeback
set yrange [] 
set ylabel '{d {/Symbol s}}/{d {p_T}} [fb]'
set format y "10^{%L}"

plot infil1 index  iindex using ($1+$2)/2:3:1:2:($3-$4):($3+$4) title mytitle

reset
set term post eps enhanced color dashed "Helvetica" 24
set output 'pwg-NLO.eps-_eta_hardest_jet_index___6.eps'
iindex = 6
mytitle = "eta_hardest_jet"
infil1 = "pwg-NLO.top"

#set log x
set log y

set datafile fortran
set style data xyerrorlines

set xtics format ""

set xrange [] writeback
set yrange [] 
set ylabel '{d {/Symbol s}}/{d {p_T}} [fb]'
set format y "10^{%L}"

plot infil1 index  iindex using ($1+$2)/2:3:1:2:($3-$4):($3+$4) title mytitle

reset
set term post eps enhanced color dashed "Helvetica" 24
set output 'pwg-NLO.eps-_pt_secondhardest_jet_index___7.eps'
iindex = 7
mytitle = "pt_secondhardest_jet"
infil1 = "pwg-NLO.top"

#set log x
set log y

set datafile fortran
set style data xyerrorlines

set xtics format ""

set xrange [] writeback
set yrange [] 
set ylabel '{d {/Symbol s}}/{d {p_T}} [fb]'
set format y "10^{%L}"

plot infil1 index  iindex using ($1+$2)/2:3:1:2:($3-$4):($3+$4) title mytitle

reset
set term post eps enhanced color dashed "Helvetica" 24
set output 'pwg-NLO.eps-_eta_secondhardest_jet_index___8.eps'
iindex = 8
mytitle = "eta_secondhardest_jet"
infil1 = "pwg-NLO.top"

#set log x
set log y

set datafile fortran
set style data xyerrorlines

set xtics format ""

set xrange [] writeback
set yrange [] 
set ylabel '{d {/Symbol s}}/{d {p_T}} [fb]'
set format y "10^{%L}"

plot infil1 index  iindex using ($1+$2)/2:3:1:2:($3-$4):($3+$4) title mytitle

