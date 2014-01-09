#!/software/bin/perl


## submitted with
## echo '~/scripts/uk10k/map_chr_pos_rsID_chrom_tabix.sh ${LSB_JOBINDEX}' | bsub -J "chr[1-22]%10" -o ~/map_chr_pos_rsID_chrom_tabix_chr%I.out


##### Collate dbSNP134 IDs for UK10K

my $chr = $ARGV[0];
# my $dbsnp="/lustre/scratch103/sanger/kw8/dbSNP/chroms/chr$chr.vcf.gz";
# my $uk10k="/lustre/scratch101/uk10k/cohorts/REL-2011-06-08/UK10K_COHORT_TWINSUK/v3/chroms/chr$chr.vcf.gz";
my $uk10k = $ARGV[1];
my $dbsnp = $ARGV[2];

## open dbSNP
open(dbSNP, "gunzip -c $dbsnp |") || die "can't open pipe to $dbsnp";
## open UK10K
open(UK10K, "gunzip -c $uk10k |") || die "can't open pipe to $uk10k";
## open outfiles
#open(OUTFILE1, ">/lustre/scratch107/projects/uk10k/cohorts/Twins/REL-2011-06-08/v3/vcf/chr_pos_rsID_2011-11-28/chr$chr".txt");
#open(OUTFILE2, ">/lustre/scratch107/projects/uk10k/cohorts/Twins/REL-2011-06-08/v3/vcf/chr_pos_rsID_2011-11-28_log/chr$chr.log");
#open(OUTFILE2, ">/lustre/scratch107/projects/uk10k/cohorts/Twins/REL-2011-06-08/v3/vcf/chr_pos_rsID_2011-12-16_log/chr$chr.log");
open(OUTFILE2, ">/lustre/scratch107/projects/uk10k/cohorts/Twins/REL-2011-06-08/v3/vcf/chr_pos_rsID_2012-01-02_log/chr$chr.log");

my $chr;
my $pos;
my $line;
my $rsid;

## open dbSNP
open(dbSNP, "gunzip -c $dbsnp |") || die "can't open pipe to $dbsnp";
## open UK10K
open(UK10K, "gunzip -c $uk10k |") || die "can't open pipe to $uk10k";

## get first entry in dbSNP
while ($rsid=<dbSNP>)  {
    ## ignore header
    if ($rsid=~/^#/)  {
    }  else  {
	chomp $rsid;
	@scols = split(" ", $rsid);
	$schr=$scols[0];
	$spos=$scols[1];
	$snpid=$scols[2];
	$sinfo=$scols[7];
        ## edit chrom (also contains MT and PAR)
	if ($schr =~ /^X/) {
	    $schr = 23;
	}
	elsif ($schr =~ /^Y/) {
	    $schr = 24;
	}
	elsif ($schr =~ /^MT/ || $schr =~ /^PAR/) {
	    $schr = 25;
	}
	print(OUTFILE2 "$schr\t$spos\tdbSNP\n");
	last;
    }
}

## loop through UK10K file
while ($line=<UK10K>)  {
    $counter=0;
    chomp $line;
    if ($line=~/^#/)  {
    }	
    else  {
	@cols = split(" ", $line);
	$chr=$cols[0];
	$pos=$cols[1];
	$info=$cols[7];
        ## edit chrom
	if ($chr =~ /^X/) {
	    $chr = 23;
	}
	elsif ($chr =~ /^Y/) {
	    $chr = 24;
	}
	print(OUTFILE2 "$chr\t$pos\tUK10K\n");
	## print rsID if coordinates match with UK10K
	if ($chr==$schr && $pos==$spos)  {
	    print("$chr\t$pos\t$snpid\t$cols[3]\t$cols[4]\t$cols[5]\t$cols[6]\t$info\t$sinfo\n");
	}
	## get next entry in dbSNP if coordinates are less than UK10K
	if ($chr==$schr && $pos>$spos || ($chr>$schr))  {
	    while ($rsid=<dbSNP>)  {
                ## ignore header
		if ($rsid=~/^#/)  {
		}
		chomp $rsid;
		@scols = split(" ", $rsid);
		$schr=$scols[0];
		$spos=$scols[1];
		$snpid=$scols[2];
		$sinfo=$scols[7];
		## edit chrom (also contains MT and PAR)
		if ($schr =~ /^X/) {
		    $schr = 23;
		}
		elsif ($schr =~ /^Y/) {
		    $schr = 24;
		}
		elsif ($schr =~ /^MT/ || $schr =~ /^PAR/) {
		    $schr = 25;
		}
		print(OUTFILE2 "$schr\t$spos\tdbSNP\n");
		## print rsID if coordinates match 
		if ($chr==$schr && $pos==$spos)  {
		    print("$chr\t$pos\t$snpid\t$cols[3]\t$cols[4]\t$cols[5]\t$cols[6]\t$info\t$sinfo\n");
		}
		## get out of loop if dbSNP coordinates are greater than UK10K
		elsif ($chr==$schr && $pos<$spos || ($chr<$schr))  {
		    $counter=1;
		    last;
		}
	    }
            ## print rsID if coordinates match with UK10K after looping through dbSNP
	    if ($chr==$schr && $pos==$spos && $counter==1)  {
		print("$chr\t$pos\t$snpid\t$cols[3]\t$cols[4]\t$cols[5]\t$cols[6]\t$info\t$sinfo\n");
	    }
	}
    }
}

close(UK10K);
close(dbSNP);
#close(OUTFILE1);
close(OUTFILE2);
