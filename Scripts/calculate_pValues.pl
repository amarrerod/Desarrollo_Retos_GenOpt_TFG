#!/usr/bin/perl
use strict;

my @init_A = qw {CMA_ES_2_50 CMA_ES_08_50 CMA_ES_03_100};
my @init_B = qw {CMA_ES_2_20 CMA_ES_08_50 CMA_ES_03_100};
my $stat;
my $pValues;

for my $problem (0..17) {
  for (my $i = 0; $i < @init_A; $i++) {
    for (my $j = $i + 1; $j < @init_B; $j++) {
      my $file_A = "/home/marrero/Escritorio/Rendimiento/CMA/2_50/boxplots/boxplots_${problem}.metric" if ($init_A[$i] eq "CMA_ES_2_50");
      $file_A = "/home/marrero/Escritorio/Rendimiento/CMA/08_50/boxplots/boxplots_${problem}.metric" if ($init_A[$i] eq "CMA_ES_08_50");
      $file_A = "/home/marrero/Escritorio/Rendimiento/CMA/03_100/boxplots/boxplots_${problem}.metric" if ($init_A[$i] eq "CMA_ES_03_100");
         

      my $file_B = "/home/marrero/Escritorio/Rendimiento/CMA/2_50/boxplots/boxplots_${problem}.metric" if ($init_A[$i] eq "CMA_ES_2_50");
      $file_B = "/home/marrero/Escritorio/Rendimiento/CMA/08_50/boxplots/boxplots_${problem}.metric" if ($init_A[$i] eq "CMA_ES_08_50");
      $file_B = "/home/marrero/Escritorio/Rendimiento/CMA/03_100/boxplots/boxplots_${problem}.metric" if ($init_A[$i] eq "CMA_ES_03_100");
             
      
      my @result = `./statisticalTests_old.pl 100 $file_A $file_B`;
      my ($mean1, $mean2, $median1, $median2, $sd1, $sd2, $pValue, $eSize) = split / /, $result[0];

      if (($i == 0) && ($j == 1)) {
        $stat->{$problem}{$init_A[$i]} = [$mean1, $median1, $sd1];
      }

      if ($i == 0) {
        $stat->{$problem}{$init_B[$j]} = [$mean2, $median2, $sd2];
      }

      if ($pValue < 0.05) {
        if (($mean1 < $mean2) && ($median1 < $median2)) {
          $pValues->{$problem}{$init_A[$i]}{$init_B[$j]} = [$pValue, '$\uparrow$'];
        }
        elsif (($mean1 > $mean2) && ($median1 > $median2)) {
          $pValues->{$problem}{$init_A[$i]}{$init_B[$j]} = [$pValue, '$\downarrow$'];
        }
      }
      else {
        $pValues->{$problem}{$init_A[$i]}{$init_B[$j]} = [$pValue, '$\leftrightarrow$'];
      }
    }
  }
}

for my $problem (0..17) {
  my $mean_min   = $stat->{$problem}{$init_A[0]}[0];
  my $median_min = $stat->{$problem}{$init_A[0]}[1];
  my $mean_alg = $init_A[0];
  my $median_alg = $init_A[0];
  my $line = "f_\{${problem}\}";

  for my $initA (@init_A) {
    if ($stat->{$problem}{$initA}[0] < $mean_min) {
      $mean_min = $stat->{$problem}{$initA}[0];
      $mean_alg = $initA;
    }
    if ($stat->{$problem}{$initA}[1] < $median_min) {
      $median_min = $stat->{$problem}{$initA}[1];
      $median_alg = $initA;
    }
  }

  for my $initA (@init_A) {
    if ($initA eq $mean_alg) {
      $line .= sprintf(' & {\\bf %.3e}', $stat->{$problem}{$initA}[0]);
    }
    else {
      $line .= sprintf(' & %.3e', $stat->{$problem}{$initA}[0]);
    }

    if ($initA eq $median_alg) {
      $line .= sprintf(' & {\\bf %.3e}', $stat->{$problem}{$initA}[1]);
    }
    else {
      $line .= sprintf(' & %.3e', $stat->{$problem}{$initA}[1]);
    }
      
    $line .= sprintf(' & %.3e', $stat->{$problem}{$initA}[2]);
  }
  $line .= ' \\\\';
  print "$line\n";
}

print "\n\n";

for my $problem (0..17) {
  my $line = "f_\{${problem}\}";
  for (my $i = 0; $i < @init_A; $i++) {
    for (my $j = $i + 1; $j < @init_B; $j++) {
      my $initA = $init_A[$i];
      my $initB = $init_B[$j];
      if ($pValues->{$problem}{$initA}{$initB}[1] ne '$\leftrightarrow$') {
       $line .= sprintf(' & {\\bf %.3e} & %s', $pValues->{$problem}{$initA}{$initB}[0], $pValues->{$problem}{$initA}{$initB}[1]);
      }
      else {
       $line .= sprintf(' & %.3e & %s', $pValues->{$problem}{$initA}{$initB}[0], $pValues->{$problem}{$initA}{$initB}[1]);
      }
    }
  }
  $line .= ' \\\\';
  print "$line\n";
}
