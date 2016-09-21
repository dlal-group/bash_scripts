#!/usr/bin/perl -w
use strict;
use Getopt::Long;


###############################################################
##################### DOCUMENTATION ###########################

## Author: Yun Li
## Date: Feb 2007
## Goal: Calculate Squared Dosage Correlation: Correlation between estimated dosages and observed genotypes
my $location="/home/ylwtx/pl/doseR2.pl";
my $title = " -- ";

print "$title\n";
print "\t\@: $location\n";
print "\tby: Yun Li(ylwtx\@umich.edu)\n\n";

## allow markers to be different (.dose and .ped)
## evaluate only those markers available in both
## match person ID's
## can handle gzipped dose and ped files
 
## Input Parameters
## 

## Input Files
## (1) dose: dosage file (2 extra columns : ID, DUMMY)
## (2) ped: observed genotypes 
## (3) dat (associated with .ped)
## (4) snps: associated with .dose
## (5) info: to read in the reference allele

## Output Files

## change log
## 1) handle indels
## 	by Qing Duan, Dec 17, 2012
## 	Note: ped allele and dose allele should both be coded in ATGC code or both in 1234 code.

################# END OF DOCUMENTATION ########################
###############################################################

my $MIN_NUM_OPTS = 6;

if($#ARGV<($MIN_NUM_OPTS*2 - 1))
{
	&usage();
	die "Failed to parse parameters\n\n";
}


my %opts = ();

	# Default Optional Options
$opts{i} = "yun.in";
$opts{o} = "yun.out";

Getopt::Long::GetOptions(\%opts,qw(
	dose=s
	ped=s
	dat=s
	snps=s
	info=s
	o=s
)) || die "Failed to parse options\n\n";


&printPar();

my %dosesnps = (); # key: marker name; # value: index (0,...);
my @doseIndices = (); # indices (0,..) in .snps for SNPs in both .dose and .ped;
my @pedIndices = (); # indices (0,..) in .dat for SNPs in both .dose and .ped;
my @snps = ();
my @refAlleles = ();
my %doseIDs = (); # key: ID; #vlaue: line index (0,..) in .dose;
my %relevantPedLines = (); # key: ID; #value: entire line;
my @sumxy = ();
my @sumx = ();
my @sumy = ();
my @sumsqx = ();
my @sumsqy = ();
my @numObsG = ();
my @doseR2 = ();

my %ACGTcodes = ( A => 1, C => 2, G => 3, T => 4, 1 => 1, 2 => 2, 3 => 3, 4 => 4);

if ($opts{dose} =~ /(.+)\.gz$/) { $opts{dose} = $1; }   # Remove gz suffix if provided
if ($opts{ped} =~ /(.+)\.gz$/) { $opts{ped} = $1; }   # Remove gz suffix if provided

&main();


# begin sub

sub usage
{
	print "\n";
	print "Usage:\n\t";
	print "-dose\t .dose file, 2 heading columns (Required) \n\t";
	print "-ped\t .ped file (Required) \n\t";
	print "-dat\t .dat file, non-M records ALLOWED (Required) \n\t";
	print "-snps\t .snps file (Required) \n\t";
	print "-info\t .info file, one heading line (Required) \n\t";
	print "-o\t Output Prefix (Required, .doseR2) \n\t";
	print "\n";
}


sub printPar
{
	print "\n";
	print "Parameters in Effect:\n";
	print "\t .dose file, 2 heading columns \n\t\t-dose '$opts{dose}'\n";
	print "\t .ped file \n\t\t-ped '$opts{ped}'\n";
	print "\t .dat file \n\t\t-dat '$opts{dat}'\n";
	print "\t .snps file \n\t\t-snps '$opts{snps}'\n";
	print "\t .info file \n\t\t-info '$opts{info}'\n";
	print "\t Output Prefix \n\t\t-o '$opts{o}'\n";
	print "\n";
}



sub rm
{
	if(-e $_[0])
	{
		print "WARNING: $_[0] existed and removed!!\n\n";
		system "/bin/rm -rf $_[0]";
	}
}



sub bakrm
{
	if(-e $_[0])
	{
		my $tt = $_[0].".old";
		print "WARNING: $_[0] existed and moved to $tt !!\n\n";
		system "mv $_[0] $tt";
	}
}


sub swap
{
	my $tmp = $_[0];
	$_[0] = $_[1];
	$_[1] = $tmp;
}


sub main
{
	&sanityChecks();
	&calcR2();
}


sub sanityChecks
{
	# (1) on markers;
	my %infosnps = ();
	open(IN, $opts{info}) || die "Can't open file $opts{info}\n\n";
	my $linenum = 0;
	while(my $line = <IN>)
	{
  	$linenum ++;
  	chomp($line);
  	$line =~ s/^\s+//;
  	my @lineArr = (split(/\s+/,$line));
  	if ($#lineArr < 1)
  	{
    	print "ERROR: .info file should have at least 2 fields, while line $linenum in $opts{info} is:\n";
    	print "       $line\n\n";
    	exit; 
  	}
  	if (!exists $infosnps{$lineArr[0]})
  	{
    	$infosnps{$lineArr[0]} = $lineArr[1];
  	}
  	else
  	{
    	print "ERROR: SNP $lineArr[0] on line $linenum in $opts{info} is duplicated\n\n";
    	exit;
  	}
	}
	close(IN);
	print "Read $linenum SNPs from $opts{info}\n\n";

	open(IN, $opts{snps}) || die "Can't open file $opts{snps}\n\n";
	$linenum = 0;
	while(my $line = <IN>)
	{
  	$linenum ++;
  	chomp($line);
  	$line =~ s/^\s+//;
  	my @lineArr = (split(/\s+/,$line));
  	if ($#lineArr != 0)
  	{
    	print "ERROR: .snps file should have exact one field, while line $linenum in $opts{snps} is:\n";
    	print "       $line\n\n";
    	exit; 
  	}
  	if (!exists $dosesnps{$line})
  	{
    	$dosesnps{$line} = $linenum - 1; 
  	}
  	else
  	{
    	print "ERROR: $line on line $linenum in $opts{snps} is duplicated\n\n";
    	exit;
  	}
	}
	close(IN);
	print "Read $linenum SNPs from $opts{snps}\n\n";
	
	open(IN, $opts{dat}) || die "Can't open file $opts{dat}\n\n";
	$linenum = 0;
	while(my $line = <IN>)
	{
  	$linenum ++;
  	chomp($line);
  	$line =~ s/^\s+//;
  	my @lineArr = (split(/\s+/,$line));
  	if ($#lineArr != 1)
  	{
    	print "ERROR: .dat file should have exact two fields, while line $linenum in $opts{dat} is:\n";
    	print "       $line\n\n";
    	exit; 
  	}
  	if (exists $dosesnps{$lineArr[1]})
  	{
    	push @doseIndices, $dosesnps{$lineArr[1]};
    	push @pedIndices, ($linenum - 1);
    	if(!exists $infosnps{$lineArr[1]})
    	{
      	print "ERROR: $lineArr[1] in .snps and .dat but not in .info\n\n";
      	exit;
    	}

## Removed ATGCcodes -- Qing Dec1712
    	push @refAlleles, $infosnps{$lineArr[1]};
    	push @snps, $lineArr[1];
  	}
	}
	close(IN);
	my $nSNPsInBoth = scalar(@doseIndices);
	print "Read $linenum SNPs from $opts{dat}, $nSNPsInBoth also in $opts{snps}\n\n";

	# check if SNPs in both .dat and .snps in .info
		
	
	# (2) on people;
	my $gzdose = $opts{dose}.".gz";
	if (-r $gzdose) 
	{
		open(IN, "gunzip -c $gzdose |") || die "Can't open compressed file $gzdose\n\n";
	}
	else
	{
		open(IN,$opts{dose}) || die "Can't open file $opts{dose}\n\n";
	}
	$linenum = 0;
	while (my $line = <IN>)
	{
  	$linenum ++;
  	chomp($line);
  	$line =~ s/^\s+//;
  	my @lineArr = (split(/\s+/,$line));
  	if ($#lineArr < 1)
  	{
    	print "ERROR: each line in $opts{dose} should have at least 2 fields, but\n";
    	print "       line $linenum does NOT\n\n";
    	exit;
  	}
  	my $doseFamID = (split('->',$lineArr[0]))[0];
		my $doseSubID = (split('->',$lineArr[0]))[1];
		my $doseID = $doseFamID."->".$doseSubID;
		if (!exists $doseIDs{$doseID})
		{
  		$doseIDs{$doseID} = $linenum - 1;
		}
		else
		{
  		print "ERROR: ID $doseID is duplicated in $opts{dose}\n\n";
  		exit;
		}
	}
	close(IN);
	print "Read $linenum people from $opts{dose}\n\n";
	
	my $gzped = $opts{ped}.".gz";
	if (-r $gzped)
	{
		open(IN, "gunzip -c $gzped |") || die "Can't open compressed file $gzped\n\n";
	}
	else
	{
		open(IN, $opts{ped}) || die "Can't open file $opts{ped}\n\n";
	}
	$linenum = 0;
	my $peopleInBoth = 0;
	while(my $line = <IN>)
	{
  	$linenum ++;
  	chomp($line);
  	$line =~ s/^\s+//;
  	my @lineArr = (split(/\s+/,$line));
  	if ($#lineArr < 4)
  	{
    	print "ERROR: each line in $opts{ped} should have at least 5 fields, but\n";
    	print "       line $linenum does NOT\n\n";
    	exit;
  	}
		my $pedID = $lineArr[0]."->".$lineArr[1];
		if (exists $doseIDs{$pedID})
		{
  		$peopleInBoth ++;
  		if (!exists $relevantPedLines{$pedID}) 
  		{
    		$relevantPedLines{$pedID} = $line;
  		}
  		else
  		{
    		if ($relevantPedLines{$pedID} ne $line)
    		{
      		print "ERROR: ID $pedID is duplicated and NOT consistent\n\n";
      		exit;
    		}
  		}
		}
		  
	}
	close(IN);
	print "Read $linenum people from $opts{ped}, $peopleInBoth also in $opts{dose}\n\n";
	
}


sub calcR2
{
	for(my $i=0; $i<=$#refAlleles; $i++)
	{
		$sumxy[$i] = 0;
		$sumx[$i] = 0;
		$sumy[$i] = 0;
		$sumsqx[$i] = 0;
		$sumsqy[$i] = 0;
		$numObsG[$i] = 0;
	}
	
        my $gzdose = $opts{dose}.".gz";
        if (-r $gzdose)
        {
                open(IN, "gunzip -c $gzdose |") || die "Can't open compressed file $gzdose\n\n";
        }
        else
        {
                open(IN,$opts{dose}) || die "Can't open file $opts{dose}\n\n";
        }
	my $linenum = 0;
	while(my $line = <IN>)
	{
  	$linenum ++;
  	chomp($line);
  	$line =~ s/^\s+//;
  	my @doseLineArr = (split(/\s+/,$line));
  	if ($#doseLineArr < 1)
  	{
    	print "ERROR: each line in $opts{dose} should have at least 2 fields, but\n";
    	print "       line $linenum does NOT\n\n";
    	exit;
  	}
  	my $doseFamID = (split('->',$doseLineArr[0]))[0];
		my $doseSubID = (split('->',$doseLineArr[0]))[1];
		my $doseID = $doseFamID."->".$doseSubID;
		if (exists $relevantPedLines{$doseID})
		{
  		my $pedLine = $relevantPedLines{$doseID};
  		my @pedLineArr = (split(/\s+/,$pedLine));
  		for(my $i=0; $i<=$#refAlleles; $i++)
  		{
  			my $pedIndex = $pedIndices[$i] + 5 ;
  			# only calculate r2 when genotype was observed;
  			if($pedLineArr[$pedIndex] ne "0/0" && $pedLineArr[$pedIndex] ne "N/N" && $pedLineArr[$pedIndex] ne "./.")
  			{
  				$numObsG[$i] ++;
  				my $doseIndex = $doseIndices[$i] + 2;
  				$sumx[$i] += $doseLineArr[$doseIndex];
  				$sumsqx[$i] += ($doseLineArr[$doseIndex] * $doseLineArr[$doseIndex]);
  			
## Removed ACGTcodes -- Qing Dec172012
  				my $pedGenoReCode = 0;
  				my @obsAlleles = (split("/",$pedLineArr[$pedIndex]));
  				if($obsAlleles[0] ne $refAlleles[$i])
  				{
  					if(exists $obsAlleles[0])
  					{
  						$pedGenoReCode ++;
  					}
  				}
  				if($obsAlleles[1] ne $refAlleles[$i])
  				{
  					if(exists $obsAlleles[1])
  					{
  						$pedGenoReCode ++;
  					}
  				
  				}
  				$sumy[$i] += $pedGenoReCode;
  				$sumsqy[$i] += ($pedGenoReCode * $pedGenoReCode);
  				
  				$sumxy[$i] += ($doseLineArr[$doseIndex] * $pedGenoReCode);
  			}
  		}
  	}
  }
	close(IN);
	
	my $of = $opts{o}.".doseR2";
	&rm($of);
	open(OUT,">>$of") || die "Can't append file $of\n\n";
	print OUT "SNP refAllele numObsGeno doseR2\n";
	for(my $i=0; $i<=$#refAlleles; $i++)
	{
		my $denomX = $numObsG[$i] * $sumsqx[$i] - ($sumx[$i] * $sumx[$i]);
		my $denomY = $numObsG[$i] * $sumsqy[$i] - ($sumy[$i] * $sumy[$i]);
		if($denomX == 0 || $denomY == 0)
		{
			print OUT "$snps[$i] $refAlleles[$i] $numObsG[$i] NA\n";
		}
		else
		{
			my $nom = $numObsG[$i] * $sumxy[$i] - ($sumx[$i] * $sumy[$i]);
			my $doseR2 = $nom * $nom;
			my $denom = $denomX * $denomY;
			$doseR2 /= $denom;
			print OUT "$snps[$i] $refAlleles[$i] $numObsG[$i] ";
			printf OUT ("%8.6f\n", $doseR2);
		}
		
	}
	close(OUT);
}



