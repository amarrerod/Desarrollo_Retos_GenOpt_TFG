#!/usr/bin/perl
use warnings;
use strict;

for my $problem (0..17) {

  my $script = qq{postscript("boxplot_CMA_ES_${problem}.eps", horizontal=FALSE, height=8, width=16, pointsize=14.9)
data1<-scan("/home/marrero/Escritorio/Rendimiento/CMA/2_50/boxplots/boxplots_${problem}.metric")
data2<-scan("/home/marrero/Escritorio/Rendimiento/CMA/2_100/boxplots/boxplots_${problem}.metric")
data3<-scan("/home/marrero/Escritorio/Rendimiento/CMA/2_20/boxplots/boxplots_${problem}.metric")
dataM<-matrix(c(data1,data2,data3), 100)
library("Rlab")
bplot(dataM, space = 0.6, labels = c("CMA_ES-2-50", "CMA_ES-2-100", "CMA_ES-2-20"), ylab = "Objective value", main = expression("Function f"[${problem}]*" - 10e6 evals. - 1 exec."))
dev.off()};

  open FILE, "> results/boxplot_CMA_ES_F${problem}.R";
  print FILE $script;
  close FILE;
}
