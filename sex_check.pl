#!/usr/bin/env perl
use strict;
use warnings;

my $bam = $ARGV[0];

chomp $bam;

my @a = split(/\//, $bam);
my $sample = $a[$#a];

my $depth_20 = 0;
open my $fh, "samtools depth -r 20 $bam | cut -f 3 |";
die "Could not open $bam\n" unless $fh;
while (<$fh>) {
	chomp;
	$depth_20 += $_;
}
close $fh;

my $depth_X = 0;
open my $fh2, "samtools depth -r X:2699521-154931043 $bam | cut -f 3 |";
die "Could not open $bam\n" unless $fh2;
while (<$fh2>) {
        chomp; 
        $depth_X += $_;
}
close $fh2;

my $depth_Y = 0;
open my $fh3, "samtools depth -r Y:2649521-59034049 $bam | cut -f 3 |";
die "Could not open $bam\n" unless $fh3;
while (<$fh3>) {
        chomp;
        $depth_Y += $_;
}
close $fh3;

my $cov_20 = $depth_20/63025520;
my $cov_X = $depth_X/(154931043-2699521);
my $cov_Y = $depth_Y/(59034049-2649521);

print "COV\t$sample\t$cov_20\t$cov_X\t$cov_Y\n";
