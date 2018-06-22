#!/usr/bin/perl
use warnings;
use strict;

for my $problem (0..17) {

  my $script = qq{set bmargin 4
#set xrange [0:362]
#set yrange [650:950]

set format y "%g"
set format x "%g"

set key top right
#set key spacing 1.5

set xlabel "Function evaluations" font "Helvetica,22"
set ylabel "Mean of the error" font "Helvetica,22"
set title  "Function f_\{${problem}\} - 10e6 evals. - 1 exec." font "Helvetica-Bold,22"
set term postscript eps enhanced color solid "Helvetica,20"

set xtics font "Helvetica,14"
#set xtics 0, 6000, 41400
set ytics font "Helvetica,20"

set output "images/mean_error_evolution_CMA_ES_F_${problem}.eps"

plot "/home/marrero/Escritorio/Rendimiento/CMA/2_50/boxplots/boxplots_${problem}.metric" with linespoints lt -1 pi -3 pt 7 ps 1.5 lc rgb "black" title "CMA_ES_2_50", "/home/marrero/Escritorio/Rendimiento/CMA/2_100/boxplots/boxplots_${problem}.metric" with linespoints lt -1 pi -3 pt 6 ps 1.5 lc rgb 'black' title "CMA_ES-2-100", "/home/marrero/Escritorio/Rendimiento/CMA/2_20/boxplots/boxplots_${problem}.metric" with linespoints lt -1 pi -3 pt 7 ps 1.5 lc rgb '#696969' title "CMA-ES_2-20"};

  open FILE, "> results/mean_error_evolution_CMA_ES_F${problem}.plot";
  print FILE $script;
  close FILE;
}
