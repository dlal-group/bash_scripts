#!/usr/bin/perl -w
use strict;
use warnings;

die "Project name not specified!\n" if (!defined $ARGV[0]);

my $input_vcf = $ARGV[0];
my $chrom  = $ARGV[1];

my $pval_cutoff;

if (defined $ARGV[2]){
	$pval_cutoff=$ARGV[2];
} else{
	$pval_cutoff="1e-2";
}

#

#$ENV{'R_LIBS'} = '/software/rarevar/plotting:$R_LIBS';


my @path = split(/\//,$input_vcf);

my $project_name = pop (@path);

if ($input_vcf =~ /\.vcf\.gz$/){
	$project_name =~ s/\.vcf\.gz$//;
} else {
	$project_name =~ s/\.vcf$//;
}

my $project = $project_name;

if (defined $chrom){
	
	if ($chrom =~ /^[1-9]$|^1[0-9]$|^2[0-3]$/){
		$project = $project.".chr$chrom";
		print "Project directory: $project/\n\n";
	} else {
		print "Chr $chrom not valid!\n";
		print "Project directories: $project*/\n\n";
	}

} else {

	print "Project directories: $project*/\n\n";

}



########  Re-submit failed jobs
print "Identifying failed jobs...\n";

my @failed_VT_dicho;
my @failed_VT_quant;
my @failed_CCRaVAT;
my @failed_Amelia;
my @failed_Ariel_dicho;
my @failed_Ariel_quant;
my @failed_aSum;
my @failed_SKAT_dicho;
my @failed_SKAT_quant;
my @failed_RR;


if (defined $chrom && $chrom =~ /^[1-9]$|^1[0-9]$|^2[0-3]$/){

	@failed_VT_dicho = `find $project/VT_Results_dicho/ -type f \\( -iname "*.out" ! -iname \".*\" \\) -exec grep -Hc Successfully {} \\; | sed s/:/\"\\t\"/ | awk '\$2!=1 {print \$1}'` if (-d "$project/VT_Results_dicho/");
	@failed_VT_quant = `find $project/VT_Results_quant/ -type f \\( -iname "*.out" ! -iname \".*\" \\) -exec grep -Hc Successfully {} \\; | sed s/:/\"\\t\"/ | awk '\$2!=1 {print \$1}'` if (-d "$project/VT_Results_quant/");
	@failed_CCRaVAT = `grep -Hc Successfully $project/CCRaVAT_Results/*.out | sed s/:/\"\\t\"/ | awk '\$2!=1 {print \$1}'` if (-d "$project/CCRaVAT_Results/");
	@failed_Amelia = `find $project/AMELIA_Results/ -type f \\( -iname "*.out" ! -iname \".*\" \\) -exec grep -Hc Successfully {} \\; | sed s/:/\"\\t\"/ | awk '\$2!=1 {print \$1}'` if (-d "$project/AMELIA_Results/");
	@failed_Ariel_dicho = `find $project/ARIEL_Results_dicho/ -type f \\( -iname "*.out" ! -iname \".*\" \\) -exec grep -Hc Successfully {} \\; | sed s/:/\"\\t\"/ | awk '\$2!=1 {print \$1}'` if (-d "$project/ARIEL_Results_dicho/");
	@failed_Ariel_quant = `find $project/ARIEL_Results_quant/ -type f \\( -iname "*.out" ! -iname \".*\" \\) -exec grep -Hc Successfully {} \\; | sed s/:/\"\\t\"/ | awk '\$2!=1 {print \$1}'` if (-d "$project/ARIEL_Results_quant/");	
	@failed_aSum = `find $project/aSum_Results/ -type f \\( -iname "*.out" ! -iname \".*\" \\) -exec grep -Hc Successfully {} \\; | sed s/:/\"\\t\"/ | awk '\$2!=1 {print \$1}'` if (-d "$project/aSum_Results/");
	@failed_SKAT_dicho = `grep -Hc Successfully $project/SKAT_Results_dicho/*.out | sed s/:/\"\\t\"/ | awk '\$2!=1 {print \$1}'` if (-d "$project/SKAT_Results_dicho/");
	@failed_SKAT_quant = `grep -Hc Successfully $project/SKAT_Results_quant/*.out | sed s/:/\"\\t\"/ | awk '\$2!=1 {print \$1}'` if (-d "$project/SKAT_Results_quant/");
	@failed_RR = `find $project/RR_Results/ -type f \\( -iname "*.out" ! -iname \".*\" \\) -exec grep -Hc Successfully {} \\; | sed s/:/\"\\t\"/ | awk '\$2!=1 {print \$1}'` if (-d "$project/RR_Results/");	

} else {

	@failed_VT_dicho = `find $project*/VT_Results_dicho/ -type f \\( -iname "*.out" ! -iname \".*\" \\) -exec grep -Hc Successfully {} \\; | sed s/:/\"\\t\"/ | awk '\$2!=1 {print \$1}'` if (glob "./$project*/VT_Results_dicho");
	@failed_VT_quant = `find $project*/VT_Results_quant/ -type f \\( -iname "*.out" ! -iname \".*\" \\) -exec grep -Hc Successfully {} \\; | sed s/:/\"\\t\"/ | awk '\$2!=1 {print \$1}'` if (glob "./$project*/VT_Results_quant");
	@failed_CCRaVAT = `grep -Hc Successfully $project*/CCRaVAT_Results/*.out | sed s/:/\"\\t\"/ | awk '\$2!=1 {print \$1}'` if (glob "./$project*/CCRaVAT_Results");
	@failed_Amelia = `find $project*/AMELIA_Results/ -type f \\( -iname "*.out" ! -iname \".*\" \\) -exec grep -Hc Successfully {} \\; | sed s/:/\"\\t\"/ | awk '\$2!=1 {print \$1}'` if (glob "./$project*/AMELIA_Results");
	@failed_Ariel_dicho = `find $project*/ARIEL_Results_dicho/ -type f \\( -iname "*.out" ! -iname \".*\" \\) -exec grep -Hc Successfully {} \\; | sed s/:/\"\\t\"/ | awk '\$2!=1 {print \$1}'` if (glob "./$project*/ARIEL_Results_dicho");
	@failed_Ariel_quant = `find $project*/ARIEL_Results_quant/ -type f \\( -iname "*.out" ! -iname \".*\" \\) -exec grep -Hc Successfully {} \\; | sed s/:/\"\\t\"/ | awk '\$2!=1 {print \$1}'` if (glob "./$project*/ARIEL_Results_quant");
	@failed_aSum = `find $project*/aSum_Results/ -type f \\( -iname "*.out" ! -iname \".*\" \\) -exec grep -Hc Successfully {} \\; | sed s/:/\"\\t\"/ | awk '\$2!=1 {print \$1}'` if (glob "./$project*/aSum_Results");
	@failed_SKAT_dicho = `grep -Hc Successfully $project*/SKAT_Results_dicho/*.out | sed s/:/\"\\t\"/ | awk '\$2!=1 {print \$1}'` if (glob "./$project*/SKAT_Results_dicho");
	@failed_SKAT_quant = `grep -Hc Successfully $project*/SKAT_Results_quant/*.out | sed s/:/\"\\t\"/ | awk '\$2!=1 {print \$1}'` if (glob "./$project*/SKAT_Results_quant");	
	@failed_RR = `find $project*/RR_Results/ -type f \\( -iname "*.out" ! -iname \".*\" \\) -exec grep -Hc Successfully {} \\; | sed s/:/\"\\t\"/ | awk '\$2!=1 {print \$1}'` if (glob "./$project*/RR_Results");
}


my @failed = (@failed_VT_dicho,@failed_VT_quant,@failed_CCRaVAT,@failed_Amelia,@failed_Ariel_dicho,@failed_Ariel_quant,@failed_aSum,@failed_SKAT_dicho,@failed_SKAT_quant,@failed_RR);

	
if (@failed){
	
	print "Jobs failed:\n @failed re-submitting...\n\n";


	foreach my $file (@failed){

		chomp($file);
				
		my @list=split(/\./,$file);
		my $index=$list[-2];
			
		open (RESUB, ">resub_failed.lsf");
		open(FAILED,$file);
		
		my $co = 0;
		
		while (<FAILED>){
			
			if ($_ =~ /LSBATCH/){
				$co = 1;
				next;
			}
			
			if ($co == 1){
				last if ($_ =~/^--/);
				if ($_ =~ /\#BSUB\s-J/){
					$_ =~ s/\[.+\]/\[$index\]/;
					print RESUB $_;
				} else {
					print RESUB $_;
				}
			}
			
		}
		
		close(FAILED);
		close (RESUB);
	
		my $hidden = $file;
		$hidden =~ s/lsf/\.lsf/;
		
		system ("mv $file $hidden");		

		system ("bsub < resub_failed.lsf");
		system ("rm resub_failed.lsf");

	}

} else {

	print "No failed job found!\n\n";
	
}


########  Read and hash regions' coordinates
print "Hashing genetic regions...\n";

my %reg_chr=();
my %reg_start=();
my %reg_end=();

use File::Find;

my $DIRLIST;

if (defined $chrom && $chrom =~ /^[1-9]$|^1[0-9]$|^2[0-3]$/){
	$DIRLIST = $project;
} else {
	$DIRLIST = ".";
}

find (\&eachFile1, $DIRLIST);

sub eachFile1 {
	
	my $filename = $_;
	my $fullpath = $File::Find::name;
	
	if ( (-e $filename) && $filename =~ /\.regions\.txt$/
		&& ($fullpath =~ /^\.\/($project)/ || $fullpath =~ /^($project)/) ){
		
		#print "$filename\t$fullpath\n";
		
		open (regions_in,$filename);					
		
		while (my $region = <regions_in>){
		
			next if ($region =~ /^\#/);
			my @row=split(/\s+/,$region);
			$row[0]=23 if ($row[0] eq "X");
			$row[0]=24 if ($row[0] eq "Y");
			$reg_chr{"$row[3]"}=$row[0];
			$reg_start{"$row[3]"}=$row[1];
			$reg_end{"$row[3]"}=$row[2];
			
		}
		
		close (regions_in);
		
	}
		
		
}


my $regions = keys %reg_start;

print "Number of regions:  $regions\n\n";


########  Retrieve p-values from individual files
print "Retrieving the tests' results...\n";

my @Amelia_results;
my @Ariel_dicho_results;
my @Ariel_quant_results;
my @CCRaVAT_results;
my @aSum_results;
my @SKAT_dicho_results;
my @SKAT_quant_results;
my @RR_results;
my @RR_labels;

if (defined $chrom && $chrom =~ /^[1-9]$|^1[0-9]$|^2[0-3]$/){

	@Amelia_results = `find $project/AMELIA_Results/ -name \"*.txt\" -exec cat {} \\;` if (-d "$project/AMELIA_Results/");
	@Ariel_dicho_results = `find $project/ARIEL_Results_dicho/ -name \"*.txt\" -exec cat {} \\;` if (-d "$project/ARIEL_Results_dicho/");
	@Ariel_quant_results = `find $project/ARIEL_Results_quant/ -name \"*.txt\" -exec cat {} \\;` if (-d "$project/ARIEL_Results_quant/");
	@CCRaVAT_results = `grep -h -vE \"Gene|Command\" $project/CCRaVAT_Results/Chr*/*_CCRVgene_MAF????.txt` if (-d "$project/CCRaVAT_Results/");
	@aSum_results = `find $project/aSum_Results/ -name \"*.txt\" -exec grep --with-filename \"\" {} \\; | sed s/:/\"\\t\"/ | sed 's/\.txt//' | sed s/\"\\/\"/\"\\t\"/g | awk '{print \$3\" \"\$4\" \"\$5\" \"\$6\" \"\$7\" \"\$8\" \"\$9\" \"\$10}'` if (-d "$project/aSum_Results/");
	@SKAT_dicho_results = `cat $project/SKAT_Results_dicho/*.results.txt` if (-d "$project/SKAT_Results_dicho/");
	@SKAT_quant_results = `cat $project/SKAT_Results_quant/*.results.txt` if (-d "$project/SKAT_Results_quant/");
	@RR_results = `find $project/RR_Results/ -name \"*.txt\" -exec cat {} \\;` if (-d "$project/RR_Results/");
	@RR_labels = `find $project/RR_Results/ -name \"*.Rout\" -exec grep -cH Error {} \\; | sed s/:/\"\\t\"/ | awk '\$2!=1' | awk 'NR==1' | awk '{print \"cut -d \\\" \\\" -f1 \"\$1}' | sh +x | grep rr` if (-d "$project/RR_Results/");
		
} else {
	
	@Amelia_results = `find $project*/AMELIA_Results/ -name \"*.txt\" -exec cat {} \\;` if (glob "./$project*/AMELIA_Results");
	@Ariel_dicho_results = `find $project*/ARIEL_Results_dicho/ -name \"*.txt\" -exec cat {} \\;` if (glob "./$project*/ARIEL_Results_dicho");
	@Ariel_quant_results = `find $project*/ARIEL_Results_quant/ -name \"*.txt\" -exec cat {} \\;` if (glob "./$project*/ARIEL_Results_quant");
	@CCRaVAT_results = `grep -h -vE \"Gene|Command\" $project*/CCRaVAT_Results/Chr*/*_CCRVgene_MAF????.txt` if (glob "./$project*/CCRaVAT_Results");
	@aSum_results = `find $project*/aSum_Results/ -name \"*.txt\" -exec grep --with-filename \"\" {} \\; | sed s/:/\"\\t\"/ | sed 's/\.txt//' | sed s/\"\\/\"/\"\\t\"/g | awk '{print \$3\" \"\$4\" \"\$5\" \"\$6\" \"\$7\" \"\$8\" \"\$9\" \"\$10}'` if (glob "./$project*/aSum_Results");
	@SKAT_dicho_results = `cat $project*/SKAT_Results_dicho/*.results.txt` if (glob "./$project*/SKAT_Results_dicho");
	@SKAT_quant_results = `cat $project*/SKAT_Results_quant/*.results.txt` if (glob "./$project*/SKAT_Results_quant");
	@RR_results = `find $project*/RR_Results/ -name \"*.txt\" -exec cat {} \\;` if (glob "./$project*/RR_Results");
	@RR_labels = `find $project*/RR_Results/ -name \"*.Rout\" -exec grep -cH Error {} \\; | sed s/:/\"\\t\"/ | awk '\$2!=1' | awk 'NR==1' | awk '{print \"cut -d \\\" \\\" -f1 \"\$1}' | sh +x | grep rr` if (glob "./$project*/RR_Results");
}


#

my @VT_dicho_results;
my @VT_quant_results;

find (\&eachFile2, $DIRLIST);

sub eachFile2 {
	
	my $filename = $_;
	my $fullpath = $File::Find::name;
	
	if ((-e $filename) && $filename =~ /\.Rout$/
			&& $fullpath =~/VT_Results_/ ){
			
		my $my_id = $filename;
		grep s/\.Rout//, $my_id;

		open (INPUT, $filename);
		
		#print "$fullpath\n";
		
		my $co 	= 0;
		
					
		my $T1 	= "-";
		my $T1P	= "-";
		my $T5 	= "-";
		my $T5P	= "-";
		my $WE 	= "-";
		my $WEP = "-";
		my $VT 	= "-";									
		my $VTP = "-";
		
		my @scores = {};
		
		while (<INPUT>){
			
		
			if ($_ =~ /p-values/){
				$co++; next;	
			}
			
			if ($co == 1){
									
				next if ($_ =~ /score/);
				
				my @row = split(/\s+/,$_);

				#print "$filename $#row $_\n";


				for (my $i = 1; $i <= $#row; $i++){
					push (@scores,$row[$i]);
				}
							
			}				
			
			
		}
		
		$T1 = $scores[1] if (defined $scores[1]);
		$T1P= $scores[2] if (defined $scores[2]);
		$T5 = $scores[3] if (defined $scores[3]);
		$T5P= $scores[4] if (defined $scores[4]);
		$WE = $scores[5] if (defined $scores[5]);
		$WEP= $scores[6] if (defined $scores[6]);
		$VT = $scores[7] if (defined $scores[7]);
		$VTP= $scores[8] if (defined $scores[8]);
		
		close (INPUT);
		
		if ($fullpath =~ /VT_Results_dicho/){
			push (@VT_dicho_results,"$my_id\t$T1\t$T1P\t$T5\t$T5P\t$WE\t$WEP\t$VT\t$VTP");
		} elsif ($fullpath =~ /VT_Results_quant/) {
			push (@VT_quant_results,"$my_id\t$T1\t$T1P\t$T5\t$T5P\t$WE\t$WEP\t$VT\t$VTP");	
		}
		
	}

}


########  print out entire p-values

my $Amelia_file;

if (defined $chrom && $chrom =~ /^[1-9]$|^1[0-9]$|^2[0-3]$/){
	$Amelia_file = "$project/Results.AMELIA.$project\.txt";
} else {
	$Amelia_file = "Results.AMELIA.$project\.txt";
}


my %KBAT_pval;
my %AMELIA_pval;

if (@Amelia_results){
	
	open (AMELIA, ">$Amelia_file");

	print AMELIA "region\tchr\tstart\tend\tn_snps\tpval_kbat\tpval_amelia\tperm\n";
	
	foreach my $result (@Amelia_results){
		
		my @row=split(/\s+/,$result);
			
		print AMELIA "$row[0]\t".$reg_chr{$row[0]}."\t".$reg_start{$row[0]}."\t".$reg_end{$row[0]}."\t"."$row[1]\t$row[2]\t$row[3]\t$row[4]\n";

		$row[2] =~ s/NA/-/;
		$row[3] =~ s/NA/-/;
	
		$KBAT_pval{$row[0]}=$row[2];
		$AMELIA_pval{$row[0]}=$row[3];
	
	}

	close (AMELIA);
	
}

#

my $Ariel_dicho_file;
my $Ariel_quant_file;

if (defined $chrom && $chrom =~ /^[1-9]$|^1[0-9]$|^2[0-3]$/){
	$Ariel_dicho_file = "$project/Results.ARIEL_dicho.$project\.txt";
	$Ariel_quant_file = "$project/Results.ARIEL_quant.$project\.txt";
} else {
	$Ariel_dicho_file = "Results.ARIEL_dicho.$project\.txt";
	$Ariel_quant_file = "Results.ARIEL_quant.$project\.txt";
}


my %ARIEL_dicho_pval;
my %ARIEL_weight_dicho_pval;

if (@Ariel_dicho_results){
	
	open (ARIEL_DICHO, ">$Ariel_dicho_file");

	print ARIEL_DICHO "region\tchr\tstart\tend\tn_snps\trare_snps\tpval_original\tOR_original\tLowerConf_o\tUpperConf_o\tpval_weighted\tOR_weighted\tLowerConf_w\tUpperConf_w\n";
	
	foreach my $result (@Ariel_dicho_results){
		
		my @row=split(/\s+/,$result);
		
		print ARIEL_DICHO "$row[0]\t".$reg_chr{$row[0]}."\t".$reg_start{$row[0]}."\t".$reg_end{$row[0]}."\t"."$row[1]\t$row[2]\t$row[3]\t$row[5]\t$row[7]\t$row[9]\t$row[4]\t$row[6]\t$row[8]\t$row[10]\n";

		$row[3] =~ s/NA/-/;
		$row[4] =~ s/NA/-/;
		
		$ARIEL_dicho_pval{$row[0]}=$row[3];
		$ARIEL_weight_dicho_pval{$row[0]}=$row[4];


	}

	close (ARIEL_DICHO);
	
}

#

my %ARIEL_quant_pval;
my %ARIEL_weight_quant_pval;

if (@Ariel_quant_results){
	
	open (ARIEL_QUANT, ">$Ariel_quant_file");

	print ARIEL_QUANT "region\tchr\tstart\tend\tn_snps\trare_snps\tpval_original\tOR_original\tLowerConf_o\tUpperConf_o\tpval_weighted\tOR_weighted\tLowerConf_w\tUpperConf_w\n";
	
	foreach my $result (@Ariel_quant_results){
		
		my @row=split(/\s+/,$result);
		
		print ARIEL_QUANT "$row[0]\t".$reg_chr{$row[0]}."\t".$reg_start{$row[0]}."\t".$reg_end{$row[0]}."\t"."$row[1]\t$row[2]\t$row[3]\t$row[5]\t$row[7]\t$row[9]\t$row[4]\t$row[6]\t$row[8]\t$row[10]\n";

		$row[3] =~ s/NA/-/;
		$row[4] =~ s/NA/-/;
		
		$ARIEL_quant_pval{$row[0]}=$row[3];
		$ARIEL_weight_quant_pval{$row[0]}=$row[4];


	}

	close (ARIEL_QUANT);
	
}


#

my $CCRaVAT_file;

if (defined $chrom && $chrom =~ /^[1-9]$|^1[0-9]$|^2[0-3]$/){
	$CCRaVAT_file = "$project/Results.CCRaVAT.$project\.txt";
} else {
	$CCRaVAT_file = "Results.CCRaVAT.$project\.txt";
}


my %CCRaVAT_Pearson;
my %CCRaVAT_FishExc;

if (@CCRaVAT_results){
	
	open (CCRaVAT, ">$CCRaVAT_file");

	print CCRaVAT "Gene\tChr\tStart\tStop\t(SNPs/MAF<cutoff)\tCase+RV\tCase-RV\tCont+RV\tCont-RV\tChisq\tPearsonPval\tFishExPval\tPermutations\n";

	foreach my $result (@CCRaVAT_results){
		
		my @row=split(/\s+/,$result);
		
		print CCRaVAT $result;
		
		$row[11]="-" if ($row[11] =~ /No/); 
		
		$CCRaVAT_Pearson{$row[0]} = $row[10] if ( $row[10]!=1 );
		$CCRaVAT_FishExc{$row[0]} = $row[11] if ( $row[10]!=1 );
		
		
	}

	close (CCRaVAT);
	
}

#

my $VT_dicho_file;
my $VT_quant_file;

if (defined $chrom && $chrom =~ /^[1-9]$|^1[0-9]$|^2[0-3]$/){
	$VT_dicho_file = "$project/Results.VT_dicho.$project\.txt";
	$VT_quant_file = "$project/Results.VT_quant.$project\.txt";
} else {
	$VT_dicho_file = "Results.VT_dicho.$project\.txt";
	$VT_quant_file = "Results.VT_quant.$project\.txt";
}



my %VTtest_dicho_FT01;
my %VTtest_dicho_FT01P;
my %VTtest_dicho_FT05;
my %VTtest_dicho_FT05P;
my %VTtest_dicho_WE;
my %VTtest_dicho_WEP;
my %VTtest_dicho_VT;
my %VTtest_dicho_VTP;

if (@VT_dicho_results){
	
	open (VT_DICHO, ">$VT_dicho_file");

	print VT_DICHO "region\tchr\tstart\tend\tpval_maf_0.01\tpval_maf_0.01_pphen\tpval_maf_0.05\tpval_maf_0.05_pphen\tpval_MB\tpval_MB_pphen\tpval_variable_maf\tpval_variable_maf_pphen\n";

	foreach my $result (@VT_dicho_results){
		
		my @row=split(/\s+/,$result);
		
		$VTtest_dicho_FT01{$row[0]}  = $row[1];
		$VTtest_dicho_FT01P{$row[0]} = $row[2];
		$VTtest_dicho_FT05{$row[0]}  = $row[3];
		$VTtest_dicho_FT05P{$row[0]} = $row[4];
		$VTtest_dicho_WE{$row[0]}    = $row[5];
		$VTtest_dicho_WEP{$row[0]}   = $row[6];
		$VTtest_dicho_VT{$row[0]}    = $row[7];
		$VTtest_dicho_VTP{$row[0]}   = $row[8];

		$VTtest_dicho_VT{$row[0]} = "-" if ($VTtest_dicho_VT{$row[0]} eq "NA");
		$VTtest_dicho_VTP{$row[0]} = "-" if ($VTtest_dicho_VTP{$row[0]} eq "NA");
		
		print VT_DICHO "$row[0]\t".$reg_chr{$row[0]}."\t".$reg_start{$row[0]}."\t".$reg_end{$row[0]}."\t"."$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]\t$row[8]\n";
	
	}

	close (VT_DICHO);
	
}


#

my %VTtest_quant_FT01;
my %VTtest_quant_FT01P;
my %VTtest_quant_FT05;
my %VTtest_quant_FT05P;
my %VTtest_quant_WE;
my %VTtest_quant_WEP;
my %VTtest_quant_VT;
my %VTtest_quant_VTP;

if (@VT_quant_results){
	
	open (VT_QUANT, ">$VT_quant_file");

	print VT_QUANT "region\tchr\tstart\tend\tpval_maf_0.01\tpval_maf_0.01_pphen\tpval_maf_0.05\tpval_maf_0.05_pphen\tpval_MB\tpval_MB_pphen\tpval_variable_maf\tpval_variable_maf_pphen\n";

	foreach my $result (@VT_quant_results){
		
		my @row=split(/\s+/,$result);
		
		$VTtest_quant_FT01{$row[0]}  = $row[1];
		$VTtest_quant_FT01P{$row[0]} = $row[2];
		$VTtest_quant_FT05{$row[0]}  = $row[3];
		$VTtest_quant_FT05P{$row[0]} = $row[4];
		$VTtest_quant_WE{$row[0]}    = $row[5];
		$VTtest_quant_WEP{$row[0]}   = $row[6];
		$VTtest_quant_VT{$row[0]}    = $row[7];
		$VTtest_quant_VTP{$row[0]}   = $row[8];

		$VTtest_quant_VT{$row[0]} = "-" if ($VTtest_quant_VT{$row[0]} eq "NA");
		$VTtest_quant_VTP{$row[0]} = "-" if ($VTtest_quant_VTP{$row[0]} eq "NA");
		
		print VT_QUANT "$row[0]\t".$reg_chr{$row[0]}."\t".$reg_start{$row[0]}."\t".$reg_end{$row[0]}."\t"."$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]\t$row[8]\n";
	
	}

	close (VT_QUANT);
	
}

#

my $aSum_file;

if (defined $chrom && $chrom =~ /^[1-9]$|^1[0-9]$|^2[0-3]$/){
	$aSum_file = "$project/Results.aSum.$project\.txt";
} else {
	$aSum_file = "Results.aSum.$project\.txt";
}

my %Score_pval;
my %SSU_pval;
my %SSUw_pval;
my %UminP_pval;
my %Sum_pval;
my %aSumP_pval;
my %aSum_pval;


if (@aSum_results){
	
	open (ASUM, ">$aSum_file");

	print ASUM "region\tchr\tstart\tend\tpval_Score_test\tpval_SSU\tpval_SSU_weighted\tpval_UminP\tpval_Sum\tpval_aSum_permuted\tpval_aSum\n";
	
	foreach my $result (@aSum_results){
		
		my @row=split(/\s+/,$result);
			
		print ASUM "$row[0]\t".$reg_chr{$row[0]}."\t".$reg_start{$row[0]}."\t".$reg_end{$row[0]}."\t"."$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]\n";
	
		$Score_pval{$row[0]}=$row[1];
		$SSU_pval{$row[0]}=$row[2];
		$SSUw_pval{$row[0]}=$row[3];
		$UminP_pval{$row[0]}=$row[4];
		$Sum_pval{$row[0]}=$row[5];
		$aSumP_pval{$row[0]}=$row[6];
		$aSum_pval{$row[0]}=$row[7];

	}

	close (ASUM);
	
}


#

my $SKAT_dicho_file;

if (defined $chrom && $chrom =~ /^[1-9]$|^1[0-9]$|^2[0-3]$/){
	$SKAT_dicho_file = "$project/Results.SKAT_dicho.$project\.txt";
} else {
	$SKAT_dicho_file = "Results.SKAT_dicho.$project\.txt";
}


my %SKAT_dicho_pval;
my %SKAT_dicho_pval_perm;

if (@SKAT_dicho_results){
	
	open (SKAT_DICHO, ">$SKAT_dicho_file");

	print SKAT_DICHO "region\tchr\tstart\tend\tSKAT_dicho_pval\tSKAT_dicho_pval_perm\n";
	
	foreach my $result (@SKAT_dicho_results){
		
		my @row=split(/\s+/,$result);
		
		print SKAT_DICHO "$row[0]\t".$reg_chr{$row[0]}."\t".$reg_start{$row[0]}."\t".$reg_end{$row[0]}."\t$row[1]\t$row[2]\n";
		
		$SKAT_dicho_pval{$row[0]}=$row[1];
		$SKAT_dicho_pval_perm{$row[0]}=$row[2];


	}

	close (SKAT_DICHO);
	
}


#

my $SKAT_quant_file;

if (defined $chrom && $chrom =~ /^[1-9]$|^1[0-9]$|^2[0-3]$/){
	$SKAT_quant_file = "$project/Results.SKAT_quant.$project\.txt";
} else {
	$SKAT_quant_file = "Results.SKAT_quant.$project\.txt";
}


my %SKAT_quant_pval;
my %SKAT_quant_pval_perm;

if (@SKAT_quant_results){
	
	open (SKAT_QUANT, ">$SKAT_quant_file");

	print SKAT_QUANT "region\tchr\tstart\tend\tSKAT_quant_pval\tSKAT_quant_pval_perm\n";
	
	foreach my $result (@SKAT_quant_results){
		
		my @row=split(/\s+/,$result);
		
		print SKAT_QUANT "$row[0]\t".$reg_chr{$row[0]}."\t".$reg_start{$row[0]}."\t".$reg_end{$row[0]}."\t$row[1]\t$row[2]\n";
		
		$SKAT_quant_pval{$row[0]}=$row[1];
		$SKAT_quant_pval_perm{$row[0]}=$row[2];


	}

	close (SKAT_QUANT);
	
}

#

my $RR_file;

if (defined $chrom && $chrom =~ /^[1-9]$|^1[0-9]$|^2[0-3]$/){
	$RR_file = "$project/Results.RR.$project\.txt";
} else {
	$RR_file = "Results.RR.$project\.txt";
}


my %RR_pvals;

my $rr_score_labels;
my $rr_pvalue_labels;

if (@RR_results){
		
	open (RR, ">$RR_file");

	chomp(@RR_labels);

	$rr_score_labels = "score_".join("\tscore_",@RR_labels);
	$rr_pvalue_labels = "pval_".join("\tpval_",@RR_labels);
	
	print RR "region\tchr\tstart\tend\t$rr_score_labels\t$rr_pvalue_labels\n";
	
	foreach my $result (@RR_results){
		
		my @row=split(/\s+/,$result);
			
		print RR "$row[0]\t".$reg_chr{$row[0]}."\t".$reg_start{$row[0]}."\t".$reg_end{$row[0]};
		
		my $pvalues = "";
		
		for (my $i = 1; $i <= $#row; $i++){
			print RR "\t".$row[$i];
			if ($i > ($#row-1)/2+1){
				$pvalues=$pvalues."$row[$i]\t";
			}
		}
		
		$RR_pvals{$row[0]}=$pvalues;
		
		print RR "\n";

	}

	close (RR);
	
}


########  print out merged files

my $Merged_dicho_file;
my $Merged_quant_file;

if (defined $chrom && $chrom =~ /^[1-9]$|^1[0-9]$|^2[0-3]$/){
	$Merged_dicho_file = "$project/Merged_Results.dicho.$project\.txt";
	$Merged_quant_file = "$project/Merged_Results.quant.$project\.txt";
} else {
	$Merged_dicho_file = "Merged_Results.dicho.$project\.txt";
	$Merged_quant_file = "Merged_Results.quant.$project\.txt";
}


my @DICHO_Results=(@Amelia_results,@Ariel_dicho_results,@CCRaVAT_results,@VT_dicho_results,@aSum_results,@SKAT_dicho_results);

if(@DICHO_Results){

open (MERGED_DICHO,">$Merged_dicho_file");

print MERGED_DICHO "region\tchr\tstart\tend\t";
	
print MERGED_DICHO "pval_KBAT\tpval_AMELIA\t" if (@Amelia_results);
print MERGED_DICHO "pval_ARIEL_original\tpval_ARIEL_weighted\t" if (@Ariel_dicho_results);
print MERGED_DICHO "pval_CCRaVAT_Pearson\tpval_CCRaVAT_FishExc\t" if (@CCRaVAT_results);
print MERGED_DICHO "pval_maf_0.01\tpval_maf_0.01_pphen\tpval_maf_0.05\tpval_maf_0.05_pphen\tpval_MB\tpval_MB_pphen\tpval_VT\tpval_VT_pphen\t" if (@VT_dicho_results);
print MERGED_DICHO "pval_Multivariate_Score_test\tpval_SSU(Sum_of_Squared_test)\tpval_SSU_weighted\tpval_UminP(Pans_Marginal_test)\tpval_Sum\tpval_aSum_permuted\tpval_aSum\t" if (@aSum_results);
print MERGED_DICHO "pval_SKAT\tpval_SKAT_perm\t" if (@SKAT_dicho_results);
print MERGED_DICHO "\n";

foreach my $region (keys %reg_start){
	
	print MERGED_DICHO "$region\t".$reg_chr{$region}."\t".$reg_start{$region}."\t".$reg_end{$region}."\t";
	
	if (@Amelia_results){
		if (exists $KBAT_pval{$region}) {print MERGED_DICHO $KBAT_pval{$region}."\t";}else {print MERGED_DICHO "-\t";};
		if (exists $AMELIA_pval{$region}) {print MERGED_DICHO $AMELIA_pval{$region}."\t";}else {print MERGED_DICHO "-\t";};
	}
	
	if (@Ariel_dicho_results){
		if (exists $ARIEL_dicho_pval{$region}) {print MERGED_DICHO $ARIEL_dicho_pval{$region}."\t";}else {print MERGED_DICHO "-\t";};
		if (exists $ARIEL_weight_dicho_pval{$region}) {print MERGED_DICHO $ARIEL_weight_dicho_pval{$region}."\t";}else {print MERGED_DICHO "-\t";};
	}
	
	if (@CCRaVAT_results){
		if (exists $CCRaVAT_Pearson{$region}) {print MERGED_DICHO $CCRaVAT_Pearson{$region}."\t";}else {print MERGED_DICHO "-\t";};
		if (exists $CCRaVAT_FishExc{$region}) {print MERGED_DICHO $CCRaVAT_FishExc{$region}."\t";}else {print MERGED_DICHO "-\t";};
	}
	
	if (@VT_dicho_results){
		if (exists $VTtest_dicho_FT01{$region}) {print MERGED_DICHO $VTtest_dicho_FT01{$region}."\t";}else {print MERGED_DICHO "-\t";};
		if (exists $VTtest_dicho_FT01P{$region}) {print MERGED_DICHO $VTtest_dicho_FT01P{$region}."\t";}else {print MERGED_DICHO "-\t";};
	
		if (exists $VTtest_dicho_FT05{$region}) {print MERGED_DICHO $VTtest_dicho_FT05{$region}."\t";}else {print MERGED_DICHO "-\t";};
		if (exists $VTtest_dicho_FT05P{$region}) {print MERGED_DICHO $VTtest_dicho_FT05P{$region}."\t";}else {print MERGED_DICHO "-\t";};

		if (exists $VTtest_dicho_WE{$region}) {print MERGED_DICHO $VTtest_dicho_WE{$region}."\t";}else {print MERGED_DICHO "-\t";};
		if (exists $VTtest_dicho_WEP{$region}) {print MERGED_DICHO $VTtest_dicho_WEP{$region}."\t";}else {print MERGED_DICHO "-\t";};

		if (exists $VTtest_dicho_VT{$region}) {print MERGED_DICHO $VTtest_dicho_VT{$region}."\t";}else {print MERGED_DICHO "-\t";};
		if (exists $VTtest_dicho_VTP{$region}) {print MERGED_DICHO $VTtest_dicho_VTP{$region}."\t";}else {print MERGED_DICHO "-\t";};
	}
	
	if (@aSum_results){
		if (exists $Score_pval{$region}) {print MERGED_DICHO $Score_pval{$region}."\t";}else {print MERGED_DICHO "-\t";};
		if (exists $SSU_pval{$region}) {print MERGED_DICHO $SSU_pval{$region}."\t";}else {print MERGED_DICHO "-\t";};
		if (exists $SSUw_pval{$region}) {print MERGED_DICHO $SSUw_pval{$region}."\t";}else {print MERGED_DICHO "-\t";};
		if (exists $UminP_pval{$region}) {print MERGED_DICHO $UminP_pval{$region}."\t";}else {print MERGED_DICHO "-\t";};
		if (exists $Sum_pval{$region}) {print MERGED_DICHO $Sum_pval{$region}."\t";}else {print MERGED_DICHO "-\t";};
		if (exists $aSumP_pval{$region}) {print MERGED_DICHO $aSumP_pval{$region}."\t";}else {print MERGED_DICHO "-\t";};
		if (exists $aSum_pval{$region}) {print MERGED_DICHO $aSum_pval{$region}."\t";}else {print MERGED_DICHO "-\t";};
	}
	
	if (@SKAT_dicho_results){
		if (exists $SKAT_dicho_pval{$region}) {print MERGED_DICHO $SKAT_dicho_pval{$region}."\t";}else {print MERGED_DICHO "-\t";};
		if (exists $SKAT_dicho_pval_perm{$region}) {print MERGED_DICHO $SKAT_dicho_pval_perm{$region}."\t";}else {print MERGED_DICHO "-\t";};
	}

	print MERGED_DICHO "\n";
}


close (MERGED_DICHO);

}

#

my @QUANT_Results = (@Ariel_quant_results,@VT_quant_results,@SKAT_quant_results,@RR_results);

if (@QUANT_Results){

open (MERGED_QUANT,">$Merged_quant_file");

print MERGED_QUANT "region\tchr\tstart\tend\t";
print MERGED_QUANT "pval_ARIEL_original\tpval_ARIEL_weighted\t" if (@Ariel_quant_results);	
print MERGED_QUANT "pval_maf_0.01\tpval_maf_0.01_pphen\tpval_maf_0.05\tpval_maf_0.05_pphen\tpval_MB\tpval_MB_pphen\tpval_VT\tpval_VT_pphen\t" if (@VT_quant_results);
print MERGED_QUANT "pval_SKAT\tpval_SKAT_perm\t" if (@SKAT_quant_results);
print MERGED_QUANT "$rr_pvalue_labels" if (@RR_results);
print MERGED_QUANT "\n";

foreach my $region (keys %reg_start){
	
	print MERGED_QUANT "$region\t".$reg_chr{$region}."\t".$reg_start{$region}."\t".$reg_end{$region}."\t";

	if (@Ariel_quant_results){
		if (exists $ARIEL_quant_pval{$region}) {print MERGED_QUANT $ARIEL_quant_pval{$region}."\t";}else {print MERGED_QUANT "-\t";};
		if (exists $ARIEL_weight_quant_pval{$region}) {print MERGED_QUANT $ARIEL_weight_quant_pval{$region}."\t";}else {print MERGED_QUANT "-\t";};
	}

	if (@VT_quant_results){
		if (exists $VTtest_quant_FT01{$region}) {print MERGED_QUANT $VTtest_quant_FT01{$region}."\t";}else {print MERGED_QUANT "-\t";};
		if (exists $VTtest_quant_FT01P{$region}) {print MERGED_QUANT $VTtest_quant_FT01P{$region}."\t";}else {print MERGED_QUANT "-\t";};
	
		if (exists $VTtest_quant_FT05{$region}) {print MERGED_QUANT $VTtest_quant_FT05{$region}."\t";}else {print MERGED_QUANT "-\t";};
		if (exists $VTtest_quant_FT05P{$region}) {print MERGED_QUANT $VTtest_quant_FT05P{$region}."\t";}else {print MERGED_QUANT "-\t";};

		if (exists $VTtest_quant_WE{$region}) {print MERGED_QUANT $VTtest_quant_WE{$region}."\t";}else {print MERGED_QUANT "-\t";};
		if (exists $VTtest_quant_WEP{$region}) {print MERGED_QUANT $VTtest_quant_WEP{$region}."\t";}else {print MERGED_QUANT "-\t";};

		if (exists $VTtest_quant_VT{$region}) {print MERGED_QUANT $VTtest_quant_VT{$region}."\t";}else {print MERGED_QUANT "-\t";};
		if (exists $VTtest_quant_VTP{$region}) {print MERGED_QUANT $VTtest_quant_VTP{$region}."\t";}else {print MERGED_QUANT "-\t";};
	}
	
	if (@SKAT_quant_results){
		if (exists $SKAT_quant_pval{$region}) {print MERGED_QUANT $SKAT_quant_pval{$region}."\t";}else {print MERGED_QUANT "-\t";};
		if (exists $SKAT_quant_pval_perm{$region}) {print MERGED_QUANT $SKAT_quant_pval_perm{$region}."\t";}else {print MERGED_QUANT "-\t";};
	}
	
	if (@RR_results){
		if (exists $RR_pvals{$region}) {
			print MERGED_QUANT $RR_pvals{$region}."\t";
		} else {
			foreach (@RR_labels){
				print MERGED_QUANT "-\t";
			}
		}
	}
	
	print MERGED_QUANT "\n";
}


close (MERGED_QUANT);

}

print "Entire p-values printed.\n\n";

######## Extract snp details of top hits

print "Printing top hits with p-value < $pval_cutoff...\n";

my $TopHits_dicho_file;
my $TopHits_quant_file;

if (defined $chrom && $chrom =~ /^[1-9]$|^1[0-9]$|^2[0-3]$/){
	$TopHits_dicho_file = "$project/TopHits.SIG$pval_cutoff.dicho.$project\.txt";
	$TopHits_quant_file = "$project/TopHits.SIG$pval_cutoff.quant.$project\.txt";
} else {
	$TopHits_dicho_file = "TopHits.SIG$pval_cutoff.dicho.$project\.txt";
	$TopHits_quant_file = "TopHits.SIG$pval_cutoff.quant.$project\.txt";
}

if (-e $Merged_dicho_file){
	system("(head -q -n1 $Merged_dicho_file; awk '{for (i=5;i<=NF;i++){if(\$i != \"-\" && \$i<$pval_cutoff){print \$0; break}} }' $Merged_dicho_file) > $TopHits_dicho_file");
}
if (-e $Merged_quant_file){
	system("(head -q -n1 $Merged_quant_file; awk '{for (i=5;i<=NF;i++){if(\$i != \"-\" && \$i<$pval_cutoff){print \$0; break}} }' $Merged_quant_file) > $TopHits_quant_file");
}

my $SINGLESNP_dicho_file="";
my $SINGLESNP_quant_file="";

if (defined $chrom && $chrom =~ /^[1-9]$|^1[0-9]$|^2[0-3]$/){
	
	if (-d "$project/SingleSNP_Results_dicho"){
		system("/software/bin/R CMD BATCH \"--args $project_name $chrom\" /software/rarevar/SingleSNP/univ_logistic_geno_quality_score_readin_circulated.r $project/SingleSNP_Results_dicho/univ_logistic_geno_quality_score_readin_circulated.Rout");
		$SINGLESNP_dicho_file = "$project/SingleSNP_Results_dicho/ordered_results_weighted_chr$chrom.txt";
	}
	
	if (-d "$project/SingleSNP_Results_quant"){
		system("/software/bin/R CMD BATCH \"--args $project_name $chrom\" /software/rarevar/SingleSNP/univ_linear_reg_geno_quality_score_readin_circulated.r $project/SingleSNP_Results_quant/univ_linear_reg_geno_quality_score_readin_circulated.Rout");
		$SINGLESNP_quant_file = "$project/SingleSNP_Results_quant/ordered_results_weighted_chr$chrom.txt";
	}

} else {
	
	if (glob "./$project*/SingleSNP_Results_dicho"){
		system("/software/bin/R CMD BATCH \"--args $project_name\" /software/rarevar/SingleSNP/univ_logistic_geno_quality_score_readin_circulated.r SingleSNP.Rout");
		$SINGLESNP_dicho_file = "SingleSNP_dicho_ordered_results_weighted.txt";
		system("(head -q -n1 $project*/SingleSNP_Results_dicho/ordered_results_weighted*.txt | sort | uniq; grep -hv SNPname $project*/SingleSNP_Results_dicho/ordered_results_weighted*.txt) > $SINGLESNP_dicho_file");
		
	}
	
	if (glob "./$project*/SingleSNP_Results_quant"){
		system("/software/bin/R CMD BATCH \"--args $project_name\" /software/rarevar/SingleSNP/univ_linear_reg_geno_quality_score_readin_circulated.r SingleSNP.Rout");	
		$SINGLESNP_quant_file = "SingleSNP_quant_ordered_results_weighted.txt";
		system("(head -q -n1 $project*/SingleSNP_Results_quant/ordered_results_weighted*.txt | sort | uniq; grep -hv SNPname $project*/SingleSNP_Results_quant/ordered_results_weighted*.txt) > $SINGLESNP_quant_file");
	}
	
}


##

my %functions;

open(FNC,"/software/rarevar/VTtest/functional_effects.txt");

while(<FNC>){
	
	my @row= split(/\s+/,$_);
	$row[0]=23 if ($row[0] eq "X");
	$row[0]=24 if ($row[0] eq "Y");
	$functions{"$row[0]:$row[1]"}=$row[2];
}

close(FNC);


##

my $TopHits_SNPs_dicho_file;
my $TopHits_SNPs_quant_file;

if (defined $chrom && $chrom =~ /^[1-9]$|^1[0-9]$|^2[0-3]$/){
	$TopHits_SNPs_dicho_file = "$project/TopHits.SNPs.SIG$pval_cutoff.dicho.$project\.txt";
	$TopHits_SNPs_quant_file = "$project/TopHits.SNPs.SIG$pval_cutoff.quant.$project\.txt";
} else {
	$TopHits_SNPs_dicho_file = "TopHits.SNPs.SIG$pval_cutoff.dicho.$project\.txt";
	$TopHits_SNPs_quant_file = "TopHits.SNPs.SIG$pval_cutoff.quant.$project\.txt";
}


my %tophits_dicho = ();

if (-f $TopHits_dicho_file){
	
	open (TOPHITS_DICHO, $TopHits_dicho_file);
	
	while (<TOPHITS_DICHO>){
		next if ($_=~/pval_/);
		my @row=split(/\s+/, $_);
		$tophits_dicho{$row[0]}=".";
	}
	
	close (TOPHITS_DICHO);
	
	#
	
	if (-e $SINGLESNP_dicho_file){
		
		open (TOPHITS_DICHO_SNPS, ">$TopHits_SNPs_dicho_file");
				
		print TOPHITS_DICHO_SNPS "region\tchr\tpos\tSNPname\trefallele\taltallele\tSNPscore\tmaf\tcounts.cases.MM\tcounts.cases.Mm\tcounts.cases.mm\tcounts.controls.MM\tcounts.controls.Mm\tcounts.controls.mm\tminorref\tbeta\tse\tchisq\tn\tOR\tLOR\tUOR\tWaldpvalue\n"; #\tfunction\n";
		
		foreach my $region (keys %tophits_dicho){	
			
			my @Gene_SNPs = `awk '\$1==$reg_chr{$region} && \$2>=$reg_start{$region} && \$2<=$reg_end{$region}' $SINGLESNP_dicho_file`;
			
			foreach my $line (@Gene_SNPs){
				
				chomp($line);
				my @row=split(/\s+/, $line);
				$line =~ s/\s/\t/g;
				
				print TOPHITS_DICHO_SNPS $region."\t".$line;
				#print TOPHITS_DICHO_SNPS "\t".$functions{"$row[0]:$row[1]"} if (exists $functions{"$row[0]:$row[1]"});
				print TOPHITS_DICHO_SNPS "\n";		
			
			}
			
#			open (SINGLESNP, $SINGLESNP_dicho_file);
#		
#			while (my $line = <SINGLESNP>){
#			
#			
#				my @row=split(/\s+/, $line);
#			
#				next if ($line =~ /refallele/);
#				next if ($row[0] ne $reg_chr{$region});
#				
#				$line =~ s/\s/\t/g;
#				
#				if 	($row[0] eq $reg_chr{$region} && 
#						$row[1]>=$reg_start{$region} &&
#							$row[1]<=$reg_end{$region}){					
#					print TOPHITS_DICHO_SNPS $region."\t".$line;
#					print TOPHITS_DICHO_SNPS $functions{"$row[0]:$row[1]"} if (exists $functions{"$row[0]:$row[1]"});
#					print TOPHITS_DICHO_SNPS "\n";		
#				}
#			
#			}
#		
#			close (SINGLESNP);
		
		}
		
		close (TOPHITS_DICHO_SNPS);
	

	}
	
}

#

my %tophits_quant = ();

if (-f $TopHits_quant_file){
	
	open (TOPHITS_QUANT, $TopHits_quant_file);
	
	while (<TOPHITS_QUANT>){
		next if ($_=~/pval_/);
		my @row=split(/\s+/, $_);
		$tophits_quant{$row[0]}=".";
	}
	
	close (TOPHITS_QUANT);
	
	#
	
	if (-e $SINGLESNP_quant_file){
		
		open (TOPHITS_QUANT_SNPS, ">$TopHits_SNPs_quant_file");
				
		print TOPHITS_QUANT_SNPS "region\tchr\tpos\tSNPname\trefallele\taltallele\tSNPscore\tmaf\tminorref\tbeta\tse\tteststat\tn\tL\tU\tpvalue\n"; #\tfunction\n";
		
		foreach my $region (keys %tophits_quant){

			my @Gene_SNPs = `awk '\$1==$reg_chr{$region} && \$2>=$reg_start{$region} && \$2<=$reg_end{$region}' $SINGLESNP_quant_file`;
			
			foreach my $line (@Gene_SNPs){
				
				chomp($line);
				my @row=split(/\s+/, $line);
				$line =~ s/\s/\t/g;
				
				print TOPHITS_QUANT_SNPS $region."\t".$line;
				#print TOPHITS_QUANT_SNPS "\t".$functions{"$row[0]:$row[1]"} if (exists $functions{"$row[0]:$row[1]"});
				print TOPHITS_QUANT_SNPS "\n";		
			
			}
						
		}
		
		close (TOPHITS_QUANT_SNPS);
	

	}
	
}



print "Single-point significance of top-hits printed.\n\n";

########  Plotting the results
print "Creating QQ and Manhattan plots...\n";

if (defined $chrom && $chrom =~ /^[1-9]$|^1[0-9]$|^2[0-3]$/){

	if (-e "Merged_Results.dicho.$project.txt"){
		system ("cd $project; /software/bin/R CMD BATCH \"--args Merged_Results.dicho.$project.txt $pval_cutoff\" /software/rarevar/plotting/MhtPlot.R");
		system ("cd $project; /software/bin/R CMD BATCH \"--args Merged_Results.dicho.$project.txt $pval_cutoff\" /software/rarevar/plotting/QQPlot.R");
	}
	
	if (-e "Merged_Results.quant.$project.txt"){
		system ("cd $project; /software/bin/R CMD BATCH \"--args Merged_Results.quant.$project.txt $pval_cutoff\" /software/rarevar/plotting/MhtPlot.R");
		system ("cd $project; /software/bin/R CMD BATCH \"--args Merged_Results.quant.$project.txt $pval_cutoff\" /software/rarevar/plotting/QQPlot.R");
	}

} else {
	
	if (-e $Merged_dicho_file){
		system ("/software/bin/R CMD BATCH \"--args $Merged_dicho_file $pval_cutoff\" /software/rarevar/plotting/MhtPlot.R");
		system ("/software/bin/R CMD BATCH \"--args $Merged_dicho_file $pval_cutoff\" /software/rarevar/plotting/QQPlot.R");	
	}
	
	if (-e $Merged_quant_file){
		system ("/software/bin/R CMD BATCH \"--args $Merged_quant_file $pval_cutoff\" /software/rarevar/plotting/MhtPlot.R");
		system ("/software/bin/R CMD BATCH \"--args $Merged_quant_file $pval_cutoff\" /software/rarevar/plotting/QQPlot.R");
	}
	
}


print "Done!\n";
