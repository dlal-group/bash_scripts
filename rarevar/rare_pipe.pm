#!/usr/bin/perl


=head1 NAME

  Pipeline for quantitative trait and case-control analysis

=head1 DESCRIPTION

  This script initiates the pipeline for rare variant analysis
  of next generation DNA sequencing data. 

  The code operates several tests of rare variant association
  including, allele matching test (AMELIA), collapsing methods
  (VT, ARIEL, SKAT, ..) and single locus tests. It takes vcf, 
  phenotype and gene region files, deals with format conversions
  and submits LSF jobs. 

  Please email comments or questions to <ym3@sanger.ac.uk>.

by Yasin Memari

=cut

use strict;
use warnings;

use Getopt::Long;
use FileHandle;

use Cwd;
my $dir = getcwd;

use lib qw(.);

use lib '/software/rarevar/lib';
#use lib '/nfs/team151/ym3/exomes/stat/Pipeline/Modules'; #unflag to test before installation

use SinglePoint;  #Modules::SinglePoint;
use VTtestFormat;
use AmeliaFormat;
use aSum;
use RRFormat;

our $INSTALL_PATH = '/software/rarevar';

my $R_PATH 	= '/software/bin/R';
my $VCFTOOLS_PATH = $INSTALL_PATH.'/bin/vcftools';
my $TABIX_PATH = $INSTALL_PATH.'/bin/tabix';
my $BGZIP_PATH = $INSTALL_PATH.'/bin/bgzip';
my $PLINK_PATH = '/software/bin/plink';

######## get command-line options
our ($input_vcf, $input_pheno, $input_regions, $tests);
our ($chrom, $cutoff, $margin, $FisherSig);
our ($perm_VT, $perm_CCRaVAT, $perm_aSum, $perm_SKAT, $perm_RR);
our ($maf_amelia, $maf_skat, $aSum_alpha0, $rr_lambda, $help);

my $args = scalar @ARGV;

GetOptions(

	'vcf=s'     	=> \$input_vcf,
	'pheno=s'  	 	=> \$input_pheno,
	'regions=s'		=> \$input_regions,
	'tests=s'		=> \$tests,
	
	'chrom=s'		=> \$chrom,
	'cutoff=s'  	=> \$cutoff,
	'margin=s'  	=> \$margin,	
	'fisher=s'  	=> \$FisherSig,	
	
	'perm-vt=s'  	=> \$perm_VT,
	'perm-ccr=s'  	=> \$perm_CCRaVAT,
	'perm-aSum=s'  	=> \$perm_aSum,
	'perm-SKAT=s'  	=> \$perm_SKAT,
	'perm-RR=s'		=> \$perm_RR,
	
	'maf-amelia=s'	=> \$maf_amelia,
	'maf-skat=s'	=> \$maf_skat,
	'alpha0=s'	    => \$aSum_alpha0,
	'rr-lambda=s'   => \$rr_lambda,
	
	'help'			=> \$help,
	
);

# set defaults
$tests	||= "single vt amelia ariel ccravat asum skat rr";
$chrom	||= "1-24";

$cutoff	||= 0.05;
$maf_amelia  ||= 1;
$maf_skat	 ||= 1;

$margin	||= 0;
$FisherSig	||= 1;

$perm_VT	   ||= 100000;
$perm_CCRaVAT  ||= 1000000;
$perm_aSum	   ||= 10000;
$perm_SKAT	   ||= 10000;
$perm_RR	   ||= 100000;

$aSum_alpha0 ||= 0.1;
$rr_lambda	 ||= "1,5,10";

$chrom = 23 if ($chrom eq "X");
$chrom = 24 if ($chrom eq "Y");

# print usage message if requested or no args supplied
if(defined($help) || !$args) {
	&usage;
	exit(0);
}


##################### Check data files consistency ####################
#######################################################################

# vcf file
my $in_vcf_handle = new FileHandle;

my @Sample_IDs;

if(defined($in_vcf_handle)) {

	die("ERROR: Could not find input file ", $input_vcf, "\n") unless -e $input_vcf;

	if ($input_vcf =~ /\.vcf\.gz$/){
		open $in_vcf_handle, "gunzip -c $input_vcf|" or die "ERROR: Could not open zipped file:\t$input_vcf\n";
	} else {
		$in_vcf_handle->open($input_vcf);
	}

}

while (<$in_vcf_handle>){
	
	chomp;

	if ($_ =~ /^\#/){
	   	if ($_ =~ /^\#\#/){
			next;
		} else {
			my @vcf_header = split (/\s+/,$_);
			@Sample_IDs = @vcf_header[9..$#vcf_header];
			last;
	   	}
	}
	
}

close ($in_vcf_handle);

# project definition

my @path = split(/\//,$input_vcf);
my $project = pop (@path);

my $vcf_suffix;

if ($input_vcf =~ /\.vcf\.gz$/){
	$project =~ s/\.vcf\.gz$//;
	$vcf_suffix = "--gzvcf";
} else {
	$project =~ s/\.vcf$//;
	$vcf_suffix = "--vcf";
}


if ($chrom eq "1-24"){
	$project = $project;
} else {$project = $project.".chr$chrom";}


our $project_dir = $project;

mkdir $project_dir or die "Cann't create the project directory: $project_dir\nPrevious data exists... remove old files first!\n";


# phenotype file
our %Pheno_quant = ();
our %Pheno_dicho = ();
our %Pheno_cases = ();
our %Pheno_controls = ();

open (PHENOTYPES,$input_pheno);

while(<PHENOTYPES>){
	
	next if ($_ =~ /^\#/);
	my @row = split(/\s+/,$_);
	$Pheno_quant{$row[0]}=$row[1] if ($row[1] !~ m/n(.*)a/i); # quantitative #&& grep {$_ eq $row[0]} @Sample_IDs)
	$row[2]="NA" unless defined $row[2];
	$Pheno_dicho{$row[0]}=$row[2] if ($row[2] !~ m/n(.*)a/i); # case-control

}

close (PHENOTYPES);

#

my $dicho = keys %Pheno_dicho;
my $quants = keys %Pheno_quant;

my $fam_file_dicho = "$project_dir/".$project.".dicho.fam";
my $fam_file_quant = "$project_dir/".$project.".quant.fam";

my $pheno_vt_dicho = "$project_dir/".$project.".pheno.vt.dicho.txt";
my $pheno_vt_quant = "$project_dir/".$project.".pheno.vt.quant.txt";

my $pheno_amelia = "$project_dir/".$project.".pheno.amelia.txt";

open (FAM_DICHO, ">$fam_file_dicho");
open (FAM_QUANT, ">$fam_file_quant");

open (PHENO_VT_DICHO, ">$pheno_vt_dicho");
open (PHENO_VT_QUANT, ">$pheno_vt_quant");

open (PHENO_AMELIA, ">$pheno_amelia");

my $co = 0;

for my $i (0 .. $#Sample_IDs){

	my $my_id = $Sample_IDs[$i];

	if (exists $Pheno_dicho{$my_id} && $Pheno_dicho{$my_id} =~ m/case/i){	

			$Pheno_cases{$my_id}=".";
			print FAM_DICHO "$my_id $my_id 0 0 0 2\n";
			print PHENO_VT_DICHO "$my_id\t1\n";
			print PHENO_AMELIA "1\n";

	} elsif (exists $Pheno_dicho{$my_id} && $Pheno_dicho{$my_id} =~ m/control/i){

	    	$Pheno_controls{$my_id}=".";
	    	print FAM_DICHO "$my_id $my_id 0 0 0 1\n";
	    	print PHENO_VT_DICHO "$my_id\t-1\n";
	    	print PHENO_AMELIA "0\n";

	} else {
			
			if ($dicho!=0){
				if ($co==0) {print "Skipping ids with unknown case/control status: "};
		    	delete($Pheno_dicho{$my_id});
		    	print $my_id.", ";
		    	$co+=1;
			}

	};


}


$co = 0;

for my $i (0 .. $#Sample_IDs){

	my $my_id = $Sample_IDs[$i];

	if (exists $Pheno_quant{$my_id}){
		
		print FAM_QUANT "$my_id $my_id 0 0 0 ".$Pheno_quant{$my_id}."\n";
		print PHENO_VT_QUANT "$my_id\t".$Pheno_quant{$my_id}."\n";
	
	} else {

		if ($quants!=0){
			if ($co==0) {print "\nSkipping ids with unknown quantitative phenotype: "};
			delete($Pheno_quant{$my_id});
			print $my_id.", ";
			$co+=1;
		}
		
	}

}

close (PHENO_AMELIA);

close (PHENO_VT_QUANT);
close (PHENO_VT_DICHO);

close (FAM_QUANT);
close (FAM_DICHO);

#

$co = 0;

my %Pheno_merged = (%Pheno_quant,%Pheno_dicho);

foreach my $my_id (keys  %Pheno_merged){

	if (!grep {$_ eq $my_id} @Sample_IDs){

		if ($co==0) {print "\nSkipping ids, Genotypes not found in vcf: "};
		delete($Pheno_quant{$my_id}) if (exists $Pheno_quant{$my_id});
		delete($Pheno_dicho{$my_id}) if (exists $Pheno_dicho{$my_id});
		delete($Pheno_cases{$my_id}) if (exists $Pheno_cases{$my_id});
		delete($Pheno_controls{$my_id}) if (exists $Pheno_controls{$my_id});
		print $my_id.", ";
		$co+=1;
		
	}
	
}


#

my $cases = keys %Pheno_cases;
my $controls = keys %Pheno_controls;
$quants = keys %Pheno_quant;

print "\n\n", $cases, " cases to consider: \n"; 
for my $key ( keys %Pheno_cases) {print "$key, ";};
print "\n", $controls, " controls to consider: \n";
for my $key ( keys %Pheno_controls) {print "$key, ";};

print "\n".$quants, " samples with quantitative phenotype\n";


if ($cases+$controls==0 && $quants==0){
die "\nNo samples left for analysis!\n";}


if ($cases+$controls!=0 && $quants!=0){

	my %Pheno_merged = (%Pheno_quant,%Pheno_dicho);
	
	foreach my $my_id (keys %Pheno_merged){
		if (!exists $Pheno_dicho{$my_id}) {
			die "\nSample with no dichotomous status: $my_id\nPlease provide both of the phenotypes for the same set of samples when running quantitative and case-control simultaneously.\n"};
		if (!exists $Pheno_quant{$my_id}) {die "\nSample with no quantitative phenotype: $my_id\nPlease provide both of the phenotypes for the same set of samples when running quantitative and case-control simultaneously.\n"};
	}
	
}



# regions file
our %region_chr = ();
our %region_start = ();
our %region_end = ();


my $project_regions = "$project_dir/".$project.".regions.txt";

open (PROJECT_REGIONS,">$project_regions");

print PROJECT_REGIONS "\#chr\tstart\tend\tregion\n";

open (REGIONS,$input_regions);

while(<REGIONS>){

	next if ($_ =~ /^\#/);

	my @row = split(/\s+/,$_);	
	
	$row[0] =~ s/chr//;
	$row[0] = 23 if ($row[0] eq "X");
	$row[0] = 24 if ($row[0] eq "Y");
	next if ($row[0] ne $chrom && $chrom ne "1-24");

	$row[1] = $row[1]-int($margin)*1000;
	$row[2] = $row[2]+int($margin)*1000;
	
	$row[1]=0 if ($row[1] < 0); 
	
	$row[0] = "X" if ($row[0] eq "23");
	$row[0] = "Y" if ($row[0] eq "24");
	$region_chr{$row[3]}=$row[0];
	$region_start{$row[3]}=$row[1];
	$region_end{$row[3]}=$row[2];

	print PROJECT_REGIONS $region_chr{$row[3]}."\t".$region_start{$row[3]}."\t".$region_end{$row[3]}."\t".$row[3]."\n";


}

close (REGIONS);
close (PROJECT_REGIONS);


my $regions_size = keys %region_start;
print "\n\n", $regions_size, " regions to consider.\n";

print "\nStarting project $project in current working directory $dir..\n\n";
print "\nRunning tests: $tests\n";

################### Extract the relevant genotypes ####################
#######################################################################
print "\nNote: Tri-allelic sites will be ignored throughout this programme.\n";

my $input_tbi;
my $regions_vcf="";

if ($input_vcf =~ /\.vcf\.gz$/){
	
	$input_tbi = $input_vcf.".tbi";
	
	if (-f $input_tbi){
		
		print "\nTabix file identified. Extracting the relevant regions...\n";
		$regions_vcf = "$project_dir/".$project.".regions.vcf.gz";
		system("(gunzip -c $input_vcf | head -1000 | grep \"^#\"; cat $project_regions | awk '{print \$1\":\"\$2\"-\"\$3}' | xargs $TABIX_PATH $input_vcf | sort -k1,1 -k2,2n | uniq) | $BGZIP_PATH > $regions_vcf");
		$input_vcf = $regions_vcf;
	
	} else {
	
		print "Tabix indexing not found! pipeline may run slower with large files.\n";
	
	}

}

#

my $ind_file;

if ($cases+$controls!=0){
	$ind_file=$fam_file_dicho;
} else {
	$ind_file=$fam_file_quant;
}

my $project_vcf = "$project_dir/".$project.".recode.vcf";

if ($chrom eq "1-24")	{
	system("$VCFTOOLS_PATH $vcf_suffix $input_vcf --keep $ind_file --bed $project_regions --min-alleles 2 --max-alleles 2 --keep-INFO CSQ --recode --out $project_dir/$project");
	system("$VCFTOOLS_PATH --vcf $project_vcf --freq --out $project_dir/$project 2>&1 1>/dev/null");
} else {
	$chrom="X" if ($chrom eq "23");
	$chrom="Y" if ($chrom eq "24");
	system("$VCFTOOLS_PATH $vcf_suffix $input_vcf --keep $ind_file --bed $project_regions --chr $chrom --min-alleles 2 --max-alleles 2 --keep-INFO CSQ --recode --out $project_dir/$project");
	system("$VCFTOOLS_PATH --vcf $project_vcf --freq --out $project_dir/$project 2>&1 1>/dev/null");
	$chrom=23 if ($chrom eq "X");
	$chrom=24 if ($chrom eq "Y");
}

if (-e $regions_vcf){
	system("rm $regions_vcf");
}

our %Minor_allele=();
our %Major_allele=();
our %Minor_freq=();
our %Major_freq=();

open (VCF_FREQ,"$project_dir/$project.frq");

while(<VCF_FREQ>){
	
	my @row=split(/\s+/,$_);
	$row[0]="23" if ($row[0] eq "X");
	$row[0]="24" if ($row[0] eq "Y");

	my @all1=split(/:/,$row[4]);

	if (!defined $row[5]){
		
		$Minor_allele{"$row[0]:$row[1]"} = ".";
		$Major_allele{"$row[0]:$row[1]"} = $all1[0];
		$Minor_freq{"$row[0]:$row[1]"} = ".";
		$Major_freq{"$row[0]:$row[1]"} = $all1[1];
	
	} else { 
	
		my @all2=split(/:/,$row[5]);	
		$Minor_allele{"$row[0]:$row[1]"} = $all1[1] <= $all2[1] ? $all1[0] : $all2[0];
		$Major_allele{"$row[0]:$row[1]"} = $all1[1] <= $all2[1] ? $all2[0] : $all1[0];
		$Minor_freq{"$row[0]:$row[1]"} = $all1[1] <= $all2[1] ? $all1[1] : $all2[1];
		$Major_freq{"$row[0]:$row[1]"} = $all1[1] <= $all2[1] ? $all2[1] : $all1[1];
	
	}

}

close (VCF_FREQ);

#### Fisher significance filtering
my %fisher_hash=();

my $fisher_tables_file;
my $fisher_pvals_file;
my $fisher_excls_file;

if ($FisherSig != 1){

	die "Option --fisher requires dichotomous status!\n" if ($cases+$controls==0);

	print "\nDetermining single point significance for exclusion\n";
	print "\nFisher's exact test (case/control) running...\n";

	$fisher_tables_file = &SinglePoint::CreateFisherTable($project_vcf);

	system("$R_PATH CMD BATCH \"--args $fisher_tables_file\" $INSTALL_PATH/SingleSNP/exactTest.R");


	#read and hash significant SNPs

	$fisher_pvals_file = $fisher_tables_file.".pvalues";
	$fisher_excls_file = $fisher_tables_file.".exclusions";

	open (FISHEROUT, ">$fisher_excls_file");

	open (FISHER, "$fisher_pvals_file");

	while (<FISHER>){
		next if ($_ =~ /chr/); 
		my @row=split(/\s+/,$_);
		if ($row[2] < $FisherSig){
			$fisher_hash{"$row[0]:$row[1]"} = $row[2];
			my $start=$row[1]-1;
			print FISHEROUT $row[0]."\t".$start."\t".$row[1]."\t".$row[2]."\n";
		}
	}

	close (FISHER);

	close (FISHEROUT);

}


#### Conversion to amelia format
my $amelia_genotypes;

our $amelia_geno_dir = "$project_dir/AMELIA_Genotypes";

if ($tests =~ m/single/i || $tests =~ m/amelia/i || $tests =~ m/ariel/i){

	$amelia_genotypes = &AmeliaFormat::FormatChange($project_vcf);

	mkdir "$amelia_geno_dir", 0777 unless -d "$amelia_geno_dir";

	foreach my $gene (keys %main::region_start){
    
    	my $chr_id = $main::region_chr{$gene};
    	
    	$chr_id = 23 if ($chr_id eq "X");
        $chr_id = 24 if ($chr_id eq "Y");

		my $geno_file = "$amelia_geno_dir/".$gene.".dat";
	
		system("awk '\$1==$chr_id && \$2>=$main::region_start{$gene} && \$2<=$main::region_end{$gene}' $amelia_genotypes > $geno_file");
	
	}

}

#### Conversion to RR format
my $rr_genotypes;

our $rr_geno_dir = "$project_dir/RR_Genotypes";

if ($tests =~ m/rr/i){

	$rr_genotypes = &RRFormat::FormatChange($project_vcf);

	mkdir "$rr_geno_dir", 0777 unless -d "$rr_geno_dir";

	foreach my $gene (keys %main::region_start){
    
    	my $chr_id = $main::region_chr{$gene};
    	
    	$chr_id = 23 if ($chr_id eq "X");
        $chr_id = 24 if ($chr_id eq "Y");

		my $geno_file = "$rr_geno_dir/".$gene.".dat";
	
		system("awk '\$1==$chr_id && \$2>=$main::region_start{$gene} && \$2<=$main::region_end{$gene}' $rr_genotypes > $geno_file");
	
	}

}


### Convert to plink formats

my $map_file = "$project_dir/$project.map";

my $dicho_files_affix = "$project_dir/$project.dicho";
my $quant_files_affix = "$project_dir/$project.quant";

my %snp_ids = ();

my $temp_files = "$project_dir/$project.temp";


if ( $tests =~ m/ccr/i || $tests =~ m/skat/i){
	
	print "\nConverting VCF to plink formats...\n";

	if ($FisherSig != 1){
		
		if ($chrom eq "1-24")	{
			system("$VCFTOOLS_PATH --vcf $project_vcf --exclude-bed $fisher_excls_file --plink --out $temp_files 2>&1 1>/dev/null");
		} else { 
			$chrom="X" if ($chrom eq "23");
			$chrom="Y" if ($chrom eq "24");
			system("$VCFTOOLS_PATH --vcf $project_vcf --chr $chrom --exclude-bed $fisher_excls_file --plink --out $temp_files 2>&1 1>/dev/null");
			$chrom=23 if ($chrom eq "X");
			$chrom=24 if ($chrom eq "Y");
		}
		
	} else {
		
		if ($chrom eq "1-24")	{
			system("$VCFTOOLS_PATH --vcf $project_vcf --plink --out $temp_files 2>&1 1>/dev/null");
		} else { 
			$chrom="X" if ($chrom eq "23");
			$chrom="Y" if ($chrom eq "24");
			system("$VCFTOOLS_PATH --vcf $project_vcf --chr $chrom --plink --out $temp_files 2>&1 1>/dev/null");
			$chrom=23 if ($chrom eq "X");
			$chrom=24 if ($chrom eq "Y");
		}	
	
	}
	
	# set the phenotypes
	
	if ($cases+$controls!=0){
		system("$PLINK_PATH --noweb --file $temp_files --pheno $fam_file_dicho --mpheno 4 --recode --out $dicho_files_affix 2>&1 1>/dev/null");
		system("mv $dicho_files_affix.map $map_file");
	}
	
	if ($quants!=0){
		system("$PLINK_PATH --noweb --file $temp_files --pheno $fam_file_quant --mpheno 4 --recode --out $quant_files_affix 2>&1 1>/dev/null");
		system("mv $quant_files_affix.map $map_file");
	}

	system("rm $temp_files.*");
	
	# fix duplicate ids;
	my $temp_map = "$temp_files.map";
	
	open (OUTMAP,">$temp_map");
	open (INMAP,"$map_file");
	
	while (<INMAP>){	
	
		my @row=split(/\s+/,$_);

		my @rsids= split(/;/,$row[1]);
		$row[1]=shift(@rsids);

		if (exists $snp_ids{"$row[1]"}){
			print OUTMAP "$row[0]\t$row[0]d$row[3]\t$row[2]\t$row[3]\n";
			#$snp_ids{"$row[0]d$row[3]"}="$row[0]:$row[3]";
		} else {
			$snp_ids{"$row[1]"}="$row[0]:$row[3]";
			print OUTMAP "$row[0]\t$row[1]\t$row[2]\t$row[3]\n";
		} 			

	}
	
	close(INMAP);
	close(OUTMAP);
	
	system("mv $temp_map $map_file");

}


################ Single-SNP linear/logistic regression ################
#######################################################################
if ($tests =~ m/single/i){

	print "\nRunning Single-point regression (original and weighted)...\n";

	our $single_result_dicho_dir = "$project_dir/SingleSNP_Results_dicho";
	our $single_result_quant_dir = "$project_dir/SingleSNP_Results_quant";

	print "Submitting jobs to the lsf..\n\n";

	my %submit_SingleSNP = (
	    bsub_opts       => "-q long -M7000000 -R 'select[mem>7000] rusage[mem=7000]'",
	    lsf_out_dicho	=> "-o $single_result_dicho_dir/lsf.%J.%I.out",
	    lsf_out_quant	=> "-o $single_result_quant_dir/lsf.%J.%I.out",
    	lsf_batch		=> "-J 'Chr$chrom-SingleSNP[$chrom]'",
    	lsf_cmd_dicho	=> "$R_PATH --vanilla --args \$LSB_JOBINDEX $amelia_genotypes $fam_file_dicho $single_result_dicho_dir/ < $INSTALL_PATH/SingleSNP/univ_logistic_geno_quality_score_circulated.r",
    	lsf_cmd_quant	=> "$R_PATH --vanilla --args \$LSB_JOBINDEX $amelia_genotypes $fam_file_quant $single_result_quant_dir/ < $INSTALL_PATH/SingleSNP/univ_linear_reg_geno_quality_score_circulated.r",
	);

	
	my $submit_single_file;

	if ($cases+$controls!=0){

		mkdir "$single_result_dicho_dir", 0777 unless -d "$single_result_dicho_dir";

		$submit_single_file = "$project_dir/submit.single.dicho.lsf";

		open (SUBMIT, ">$submit_single_file");

		print SUBMIT "\#!/bin/csh\n";
		print SUBMIT "\#BSUB $submit_SingleSNP{bsub_opts}\n";
		print SUBMIT "\#BSUB $submit_SingleSNP{lsf_out_dicho}\n";
		print SUBMIT "\#BSUB $submit_SingleSNP{lsf_batch}\n\n";

		print SUBMIT "$submit_SingleSNP{lsf_cmd_dicho}\n";

		close (SUBMIT);

		system ("bsub < $submit_single_file");

	}
	
	if ($quants!=0){
		
		mkdir "$single_result_quant_dir", 0777 unless -d "$single_result_quant_dir";

		$submit_single_file = "$project_dir/submit.single.quant.lsf";

		open (SUBMIT, ">$submit_single_file");

		print SUBMIT "\#!/bin/csh\n";
		print SUBMIT "\#BSUB $submit_SingleSNP{bsub_opts}\n";
		print SUBMIT "\#BSUB $submit_SingleSNP{lsf_out_quant}\n";
		print SUBMIT "\#BSUB $submit_SingleSNP{lsf_batch}\n\n";

		print SUBMIT "$submit_SingleSNP{lsf_cmd_quant}\n";

		close (SUBMIT);

		system ("bsub < $submit_single_file");			
			
	}
	
}

###################### Variable threshold test ########################
#######################################################################
if ($tests =~ m/vt/i){

	print "\nRunning the collpasing Variable Threshold (VT) test..\n\n";

	our $VT_geno_dir = "$project_dir/VT_AlleleCounts";
	our $VT_weight_dir = "$project_dir/VT_Weights";

	#make directories unless they already exist
	mkdir "$VT_geno_dir", 0777 unless -d "$VT_geno_dir";
	mkdir "$VT_weight_dir", 0777 unless -d "$VT_weight_dir";

	##&VTtest::vcf2cgeno($project_vcf);
	&VTtestFormat::AlleleCounts($project_vcf);

	print "Directories for VT tests created: $VT_geno_dir, $VT_weight_dir\n";
	
	#
	
	our $VT_result_dicho_dir = "$project_dir/VT_Results_dicho";
	our $VT_result_quant_dir = "$project_dir/VT_Results_quant";
	
	print "Submitting job array to the lsf queue..\n\n";

	my %submit_vt = (
    	bsub_opts       => "-q normal -M2000000 -R 'select[mem>2000] rusage[mem=2000]'",
    	lsf_out_dicho	=> "-o $VT_result_dicho_dir/lsf.%J.%I.out",
    	lsf_out_quant	=> "-o $VT_result_quant_dir/lsf.%J.%I.out",
    	lsf_batch		=> "-J 'Chr$chrom-VT[1-$regions_size]%500'",
    	lsf_cmd_dicho	=> "perl $INSTALL_PATH/VTtest/start_VT.pl \$LSB_JOBINDEX $perm_VT $pheno_vt_dicho $VT_weight_dir $VT_geno_dir $VT_result_dicho_dir",
    	lsf_cmd_quant	=> "perl $INSTALL_PATH/VTtest/start_VT.pl \$LSB_JOBINDEX $perm_VT $pheno_vt_quant $VT_weight_dir $VT_geno_dir $VT_result_quant_dir"
	);


	my $submit_VT_file;

	if ($cases+$controls!=0){
		
		mkdir "$VT_result_dicho_dir", 0777 unless -d "$VT_result_dicho_dir";

		$submit_VT_file = "$project_dir/submit.vt.dicho.lsf";

		open (SUBMIT, ">$submit_VT_file");

		print SUBMIT "\#!/bin/csh\n";
		print SUBMIT "\#BSUB $submit_vt{bsub_opts}\n";
		print SUBMIT "\#BSUB $submit_vt{lsf_out_dicho}\n";
		print SUBMIT "\#BSUB $submit_vt{lsf_batch}\n\n";

		print SUBMIT "$submit_vt{lsf_cmd_dicho}\n";

		close (SUBMIT);

		system ("bsub < $submit_VT_file");
		
	}


	if ($quants!=0){
		
		mkdir "$VT_result_quant_dir", 0777 unless -d "$VT_result_quant_dir";

		$submit_VT_file = "$project_dir/submit.vt.quant.lsf";

		open (SUBMIT, ">$submit_VT_file");

		print SUBMIT "\#!/bin/csh\n";
		print SUBMIT "\#BSUB $submit_vt{bsub_opts}\n";
		print SUBMIT "\#BSUB $submit_vt{lsf_out_quant}\n";
		print SUBMIT "\#BSUB $submit_vt{lsf_batch}\n\n";

		print SUBMIT "$submit_vt{lsf_cmd_quant}\n";

		close (SUBMIT);

		system ("bsub < $submit_VT_file");
		
	}	

} 


######### CCRaVAT (Case-Control Rare Variant Analysis Tool) ###########
#######################################################################

if ($tests =~ m/ccr/i && $cases+$controls!=0){
	
	print "\nRunning Case-Control Rare Variant Analysis Tool (CCRaVAT)..\n\n";

	our $CCRaVAT_result_dir = "$project_dir/CCRaVAT_Results";
	mkdir "$CCRaVAT_result_dir", 0777 unless -d "$CCRaVAT_result_dir";


	print "\nCreating input files for CCRaVAT...\n";

	system("perl $INSTALL_PATH/CCRaVAT/Split_GW_PedMapFiles-0.1.pl -p $dicho_files_affix.ped $map_file $project_regions");

	if ($chrom eq "1-24") {
		system("mv Chr* $CCRaVAT_result_dir");
	} else {
		if ($chrom < 10) {
			system("mv Chr0$chrom $CCRaVAT_result_dir");
		} else {
			system("mv Chr$chrom $CCRaVAT_result_dir");	
		}
	}

	##

	print "\nSubmitting job array to the lsf queue..\n";

	my %submit_CCRaVAT = (
	    bsub_opts       => "-q long -M4000000 -R 'select[mem>4000] rusage[mem=4000]'",
	    lsf_out			=> "-o $CCRaVAT_result_dir/lsf.%J.%I.out",
	    lsf_batch		=> "-J 'CCRaVAT[$chrom]'",
	    lsf_cmd			=> "perl $INSTALL_PATH/CCRaVAT/CCRaVAT-0.2.pl -gene -nchr -cstart \$LSB_JOBINDEX -cstop \$LSB_JOBINDEX -maf=$cutoff -ext=0 -pout=0.01 -fout=0.01 -pperm=0.01 -nperm=$perm_CCRaVAT -qperm -graph"
	);

	my $submit_CCRaVAT_file = "$project_dir/submit.ccravat.lsf";

	open (SUBMIT, ">$submit_CCRaVAT_file");

	print SUBMIT "\#!/bin/csh\n";
	print SUBMIT "\#BSUB $submit_CCRaVAT{bsub_opts}\n";
	print SUBMIT "\#BSUB $submit_CCRaVAT{lsf_out}\n";
	print SUBMIT "\#BSUB $submit_CCRaVAT{lsf_batch}\n\n";

	print SUBMIT "cd $CCRaVAT_result_dir/\n";
	print SUBMIT "$submit_CCRaVAT{lsf_cmd}\n";

	close (SUBMIT);

	system ("bsub < $submit_CCRaVAT_file");

} else {
	
	if ($cases+$controls==0){
		print "\nSkipping CCRaVAT, no cases/controls\n";
	};
	
}



############# AMELIA (Allele matching case-control tests) #############
#######################################################################

if ($tests =~ m/amelia/i && $cases+$controls!=0){

	print "\nAllele matching AMELIA/KBAT tests running..\n";

	our $amelia_results_dir = "$project_dir/AMELIA_Results";
	mkdir "$amelia_results_dir", 0777 unless -d "$amelia_results_dir";


	print "Result directory: $amelia_results_dir\n" ;
	print "Submitting jobs to the lsf..\n\n";

	my %submit_amelia = (
	    bsub_opts       => "-q long -M2000000 -R 'select[mem>2000] rusage[mem=2000]'",
	    lsf_out			=> "-o $amelia_results_dir/lsf.%J.%I.out",
	    lsf_batch		=> "-J 'Chr$chrom-amelia[1-$regions_size]%500'",
	    #lsf_cmd			=> "$R_PATH CMD BATCH --vanilla --no-save --max-ppsize=500000 \"--args \$LSB_JOBINDEX $maf_amelia ../../$amelia_genotypes ../../$fam_file_dicho ../../$project_regions $project.p\$LSB_JOBINDEX.txt\" $INSTALL_PATH/AMELIA/amelia_script.R $project.p\$LSB_JOBINDEX.Rout"
		lsf_cmd			=> "perl $INSTALL_PATH/AMELIA/start_AMELIA.pl \$LSB_JOBINDEX $maf_amelia $amelia_geno_dir $fam_file_dicho $project_regions $amelia_results_dir"
	);

	my $submit_amelia_file = "$project_dir/submit.amelia.lsf";

	open (SUBMIT, ">$submit_amelia_file");

	print SUBMIT "\#!/bin/csh\n";
	print SUBMIT "\#BSUB $submit_amelia{bsub_opts}\n";
	print SUBMIT "\#BSUB $submit_amelia{lsf_out}\n";
	print SUBMIT "\#BSUB $submit_amelia{lsf_batch}\n\n";

	#print SUBMIT "cd $amelia_results_dir/\n";
	print SUBMIT "$submit_amelia{lsf_cmd}\n";

	close (SUBMIT);

	system ("bsub < $submit_amelia_file");

} else {
	
	if ($cases+$controls==0){
		print "\nSkipping AMELIA, no cases/controls\n";
	};
	
}




############# ARIEL (Collapsing, minor allele proportion) #############
#######################################################################

if ($tests =~ m/ariel/i){

	print "\nCollapsing ARIEL test running..\n";

	our $ariel_result_dicho_dir = "$project_dir/ARIEL_Results_dicho";
	our $ariel_result_quant_dir = "$project_dir/ARIEL_Results_quant";

	print "Submitting jobs to the lsf..\n\n";

	my %submit_ariel = (
	    bsub_opts       => "-q normal -M2000000 -R 'select[mem>2000] rusage[mem=2000]'",
    	lsf_out_dicho	=> "-o $ariel_result_dicho_dir/lsf.%J.%I.out",
    	lsf_out_quant	=> "-o $ariel_result_quant_dir/lsf.%J.%I.out",
    	lsf_batch		=> "-J 'Chr$chrom-ariel[1-$regions_size]%500'",
    	#lsf_cmd			=> "$R_PATH CMD BATCH --vanilla --no-save --max-ppsize=500000 \"--args \$LSB_JOBINDEX $cutoff ../../$amelia_genotypes ../../$fam_file_dicho ../../$project_regions $project.p\$LSB_JOBINDEX.txt\" $INSTALL_PATH/ARIEL/ARIELqcc_script_general.R $project.chr\$LSB_JOBINDEX.Rout"
    	lsf_cmd_dicho	=> "perl $INSTALL_PATH/ARIEL/start_ARIEL.pl \$LSB_JOBINDEX $cutoff $amelia_geno_dir $fam_file_dicho $project_regions $ariel_result_dicho_dir",
 	    lsf_cmd_quant	=> "perl $INSTALL_PATH/ARIEL/start_ARIELqQT.pl \$LSB_JOBINDEX $cutoff $amelia_geno_dir $fam_file_quant $project_regions $ariel_result_quant_dir"
	);


	my $submit_ariel_file;

	if ($cases+$controls!=0){
		
		mkdir "$ariel_result_dicho_dir", 0777 unless -d "$ariel_result_dicho_dir";

		$submit_ariel_file = "$project_dir/submit.ariel.dicho.lsf";

		open (SUBMIT, ">$submit_ariel_file");

		print SUBMIT "\#!/bin/csh\n";
		print SUBMIT "\#BSUB $submit_ariel{bsub_opts}\n";
		print SUBMIT "\#BSUB $submit_ariel{lsf_out_dicho}\n";
		print SUBMIT "\#BSUB $submit_ariel{lsf_batch}\n\n";

		print SUBMIT "$submit_ariel{lsf_cmd_dicho}\n";

		close (SUBMIT);

		system ("bsub < $submit_ariel_file");
		
	}


	if ($quants!=0){
		
		mkdir "$ariel_result_quant_dir", 0777 unless -d "$ariel_result_quant_dir";

		$submit_ariel_file = "$project_dir/submit.ariel.quant.lsf";

		open (SUBMIT, ">$submit_ariel_file");

		print SUBMIT "\#!/bin/csh\n";
		print SUBMIT "\#BSUB $submit_ariel{bsub_opts}\n";
		print SUBMIT "\#BSUB $submit_ariel{lsf_out_quant}\n";
		print SUBMIT "\#BSUB $submit_ariel{lsf_batch}\n\n";

		print SUBMIT "$submit_ariel{lsf_cmd_quant}\n";

		close (SUBMIT);

		system ("bsub < $submit_ariel_file");
		
	}	


	
}


############## SNP-set Kernel Association Test (SKAT) #################
#######################################################################


if ($tests =~ m/skat/i){

	print "\nRunning SNP-set Kernel Association Test (SKAT)..\n\n";
	
	#our $SetID_file = "$project_dir/$project.SetID";
	our $SKAT_SetIDs_dir = "$project_dir/SKAT_SetIDs";

	if ($chrom eq "1-24"){

		mkdir "$SKAT_SetIDs_dir", 0777 unless -d "$SKAT_SetIDs_dir";
		
		for my $chr (1 .. 24){
			
			my $SetID_file = "$SKAT_SetIDs_dir/$project.chr$chr.SetID";
			
			open (NewSetIDFile1, ">$SetID_file");  # Creates an empty file
			
			foreach my $reg (keys %region_chr){
				next if ($region_chr{$reg} != $chr);
				system("awk '\$1==$region_chr{$reg} && \$4>=$region_start{$reg} && \$4<=$region_end{$reg} {print \"$reg \"\$2}' $map_file >> $SetID_file");
			
			}
		
		}
				
	} else {
		
		my $SetID_file = "$project_dir/$project.SetID";
		
		open (NewSetIDFile2, ">$SetID_file");  # Creates an empty file
		
		foreach my $reg (keys %region_chr){
	    	#my $SetID_file = "$SKAT_SetIDs_dir/$reg.SetID";
	    	#system("awk '\$4>=$region_start{$reg} && \$4<=$region_end{$reg} {print \"$reg \"\$2}' $map_file > $SetID_file");
	    
			system("awk '\$4>=$region_start{$reg} && \$4<=$region_end{$reg} {print \"$reg \"\$2}' $map_file >> $SetID_file");
		}
		
	}	
		
	#
	
	print "Submitting job array to the lsf queue..\n\n";

	our $SKAT_result_dicho_dir = "$project_dir/SKAT_Results_dicho";
	our $SKAT_result_quant_dir = "$project_dir/SKAT_Results_quant";

	my %submit_skat = (
    	bsub_opts       => "-q long -M7000000 -R 'select[mem>7000] rusage[mem=7000]'",
    	lsf_out_dicho	=> "-o $SKAT_result_dicho_dir/lsf.%J.%I.out",
    	lsf_out_quant	=> "-o $SKAT_result_quant_dir/lsf.%J.%I.out",
    	lsf_batch		=> "-J 'Chr$chrom-SKAT[$chrom]'", #[1-$regions_size]%500'",
    	#lsf_cmd_dicho	=> "perl $INSTALL_PATH/SKAT/start_SKAT.pl \$LSB_JOBINDEX $dicho_files_affix $fam_file_dicho $SKAT_SetIDs_dir $SKAT_result_dicho_dir",
    	#lsf_cmd_quant	=> "perl $INSTALL_PATH/SKAT/start_SKAT.pl \$LSB_JOBINDEX $quant_files_affix $fam_file_quant $SKAT_SetIDs_dir $SKAT_result_quant_dir"
    	lsf_cmd_dicho	=> "$R_PATH CMD BATCH --vanilla --no-save --max-ppsize=500000 \"--args D $perm_SKAT ../../$dicho_files_affix.bed ../../$dicho_files_affix.bim ../../$fam_file_dicho ../$project.SetID $project.dicho.SSD $project.dicho.SSD.info\" $INSTALL_PATH/SKAT/run_SKAT_Resampling.R $project.Rout",
    	lsf_cmd_quant	=> "$R_PATH CMD BATCH --vanilla --no-save --max-ppsize=500000 \"--args C $perm_SKAT ../../$quant_files_affix.bed ../../$quant_files_affix.bim ../../$fam_file_quant ../$project.SetID $project.quant.SSD $project.quant.SSD.info\" $INSTALL_PATH/SKAT/run_SKAT_Resampling.R $project.Rout",
    	lsf_cmd_dicho_nochr	=> "$R_PATH CMD BATCH --vanilla --no-save --max-ppsize=500000 \"--args D $perm_SKAT ../../$dicho_files_affix.bed ../../$dicho_files_affix.bim ../../$fam_file_dicho ../SKAT_SetIDs/$project.chr\$LSB_JOBINDEX.SetID $project.chr\$LSB_JOBINDEX.dicho.SSD $project.chr\$LSB_JOBINDEX.dicho.SSD.info\" $INSTALL_PATH/SKAT/run_SKAT_Resampling.R $project.chr\$LSB_JOBINDEX.Rout",
    	lsf_cmd_quant_nochr	=> "$R_PATH CMD BATCH --vanilla --no-save --max-ppsize=500000 \"--args C $perm_SKAT ../../$quant_files_affix.bed ../../$quant_files_affix.bim ../../$fam_file_quant ../SKAT_SetIDs/$project.chr\$LSB_JOBINDEX.SetID $project.chr\$LSB_JOBINDEX.quant.SSD $project.chr\$LSB_JOBINDEX.quant.SSD.info\" $INSTALL_PATH/SKAT/run_SKAT_Resampling.R $project.chr\$LSB_JOBINDEX.Rout"
	);
	
	#
	
	my $submit_SKAT_file;

	if ($cases+$controls!=0){

		system("$PLINK_PATH --noweb --ped $dicho_files_affix.ped --map $map_file --max-maf $maf_skat --make-bed --out $dicho_files_affix 2>&1 1>/dev/null");
		
		mkdir "$SKAT_result_dicho_dir", 0777 unless -d "$SKAT_result_dicho_dir";

		$submit_SKAT_file = "$project_dir/submit.skat.dicho.lsf";

		open (SUBMIT, ">$submit_SKAT_file");

		print SUBMIT "\#!/bin/csh\n";
		print SUBMIT "\#BSUB $submit_skat{bsub_opts}\n";
		print SUBMIT "\#BSUB $submit_skat{lsf_out_dicho}\n";
		print SUBMIT "\#BSUB $submit_skat{lsf_batch}\n\n";

		print SUBMIT "cd $SKAT_result_dicho_dir/\n";
		
		if ($chrom eq "1-24"){
			print SUBMIT "$submit_skat{lsf_cmd_dicho_nochr}\n";
		} else {
			print SUBMIT "$submit_skat{lsf_cmd_dicho}\n";
		}
		
		close (SUBMIT);

		system ("bsub < $submit_SKAT_file");
		
	}


	if ($quants!=0){

		system("$PLINK_PATH --noweb --ped $quant_files_affix.ped --map $map_file --max-maf $maf_skat --make-bed --out $quant_files_affix 2>&1 1>/dev/null");
				
		mkdir "$SKAT_result_quant_dir", 0777 unless -d "$SKAT_result_quant_dir";

		$submit_SKAT_file = "$project_dir/submit.skat.quant.lsf";

		open (SUBMIT, ">$submit_SKAT_file");

		print SUBMIT "\#!/bin/csh\n";
		print SUBMIT "\#BSUB $submit_skat{bsub_opts}\n";
		print SUBMIT "\#BSUB $submit_skat{lsf_out_quant}\n";
		print SUBMIT "\#BSUB $submit_skat{lsf_batch}\n\n";

		print SUBMIT "cd $SKAT_result_quant_dir/\n";
		
		if ($chrom eq "1-24"){
			print SUBMIT "$submit_skat{lsf_cmd_quant_nochr}\n";
		} else {
			print SUBMIT "$submit_skat{lsf_cmd_quant}\n";
		}

		close (SUBMIT);

		system ("bsub < $submit_SKAT_file");
		
	}	
	
} 


############## Adaptive Sum Test (Han & Pan) #################
#######################################################################

if ($tests =~ m/aSum/i && $cases+$controls!=0){

	print "\nRunning the Adaptive Sum Test..\n\n";

	our $aSum_geno_dir = "$project_dir/aSum_Genotypes";
	#make directories unless they already exist
	mkdir "$aSum_geno_dir", 0777 unless -d "$aSum_geno_dir";

	&aSum::vcf2geno($project_vcf);

	our $aSum_results_dir = "$project_dir/aSum_Results";
	mkdir "$aSum_results_dir", 0777 unless -d "$aSum_results_dir";

	print "Allele counts for aSum test generated: $aSum_geno_dir\n";

	print "Submitting jobs to the lsf..\n\n";

	my %submit_aSum = (
	    bsub_opts       => "-q normal -M2000000 -R 'select[mem>2000] rusage[mem=2000]'",
    	lsf_out			=> "-o $aSum_results_dir/lsf.%J.%I.out",
    	lsf_batch		=> "-J 'Chr$chrom-aSum[1-$regions_size]%500'",
    	lsf_cmd			=> "perl $INSTALL_PATH/aSum/start_aSum.pl \$LSB_JOBINDEX $pheno_amelia $aSum_geno_dir $aSum_results_dir $aSum_alpha0 $perm_aSum"
	);

	my $submit_aSum_file = "$project_dir/submit.aSum.lsf";

	open (SUBMIT, ">$submit_aSum_file");

	print SUBMIT "\#!/bin/csh\n";
	print SUBMIT "\#BSUB $submit_aSum{bsub_opts}\n";
	print SUBMIT "\#BSUB $submit_aSum{lsf_out}\n";
	print SUBMIT "\#BSUB $submit_aSum{lsf_batch}\n\n";

	print SUBMIT "$submit_aSum{lsf_cmd}\n";

	close (SUBMIT);

	system ("bsub < $submit_aSum_file");


} else {
	
	if ($cases+$controls==0){
		print "\nSkipping aSum test, no cases/controls\n";
	};
	
}


################# RR Test (Ridge Regression analysis) #################
#######################################################################

if ($tests =~ m/rr/i && $quants!=0){

	print "\nRidge Regression running..\n";

	our $rr_results_dir = "$project_dir/RR_Results";
	mkdir "$rr_results_dir", 0777 unless -d "$rr_results_dir";


	print "Result directory: $rr_results_dir\n" ;
	print "Submitting jobs to the lsf..\n\n";

	my %submit_rr = (
	    bsub_opts       => "-q normal -M2000000 -R 'select[mem>2000] rusage[mem=2000]'",
	    lsf_out			=> "-o $rr_results_dir/lsf.%J.%I.out",
	    lsf_batch		=> "-J 'Chr$chrom-RR[1-$regions_size]%500'",
		lsf_cmd			=> "perl $INSTALL_PATH/RR/start_RR.pl \$LSB_JOBINDEX $rr_geno_dir $fam_file_quant $rr_lambda $project_regions $rr_results_dir $perm_RR"
	);

	my $submit_rr_file = "$project_dir/submit.rr.lsf";

	open (SUBMIT, ">$submit_rr_file");

	print SUBMIT "\#!/bin/csh\n";
	print SUBMIT "\#BSUB $submit_rr{bsub_opts}\n";
	print SUBMIT "\#BSUB $submit_rr{lsf_out}\n";
	print SUBMIT "\#BSUB $submit_rr{lsf_batch}\n\n";

	#print SUBMIT "cd $rr_results_dir/\n";
	print SUBMIT "$submit_rr{lsf_cmd}\n";

	close (SUBMIT);

	system ("bsub < $submit_rr_file");

} else {
	
	if ($quants==0){
		print "\nSkipping Ridge Regression, no quantitative pheno.\n";
	};
	
}



######## set usage guide
sub usage {
	my $usage =<<END;
	
Usage:
perl rare_pipe.pl [arguments]

Options
-h	| --help	Display this message and quit

-v	| --vcf		Input vcf file [should be in gzipped/unzipped VCF format v4.0]
-ph	| --pheno	Input phenotype file [format: id  quant_pheno  case/control]
-r	| --regions	Input regions file [in bed format, without header]

-t	| --tests	Methods to run including:   [default runs all the tests]
			SingleSNP, AMELIA, ARIEL, VT, CCRaVAT, aSum, SKAT, RR

-chr	| --chrom	Chromosome for analysis [default: all chromosomes]
-mar	| --margin	Flanking margin to be considered, in kb [default : "0"]
-sig	| --fisher	Maximum Fisher's exact test p-value for SNP filtering (case-control) [default: "1"]

-maf	| --cutoff	Maximum maf cutoff for the collapsing tests ARIEL, CCRaVAT [default: "0.05"]
-mafam	| --maf-amelia	Maximum maf cutoff for the allele matching test AMELIA [default: "1"]
-mafsk	| --maf-skat	Maximum maf cutoff for SKAT [default: "1"]

-pvt	| --perm-vt	Number of permutations used in VT test [default: "100,000"]
-pccr	| --perm-ccr	Number of permutations used in CCRaVAT test [default: "1,000,000"]
-pasum	| --perm-aSum	Number of permutations used in aSum test [default: "10,000"]
-pskat	| --perm-SKAT	Number of re-samplings for the SKAT test [default: "10,000"]
-prr	| --perm-RR	Number of re-samplings for the RR test [default: "100,000"]

-alph	| --alpha0	Cut-off used to flip the coding in the aSum test [default: "0.1"]
-lamb	| --rr-lambda	Lambda parameters for ridge regression (RR) test [default: "1,5,10"]

END

	print $usage;
}
