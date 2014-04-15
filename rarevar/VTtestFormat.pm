package VTtestFormat;

sub AlleleCounts {

	my $input = shift; #input file in VCF format

	my $vt_format_file = $input; #pop (@path);
	$vt_format_file =~ s/\.recode\.vcf/\.VT_format\.txt/;

	open(OUTPUT,">$vt_format_file") or die "can't open output: $vt_format_file\n";

	# read polyphen weighted probs

	my %Weight_prob;
	
	my $INSTALL_PATH = $main::INSTALL_PATH;

	open (PolyPh, "$INSTALL_PATH/VTtest/VT_weights.txt");

	while (<PolyPh>){
		my @row=split(/\s+/,$_);
		my @pos=split(/_/,$row[0]);
		$pos[0]=23 if ($pos[0] eq "X");
		$pos[0]=24 if ($pos[0] eq "Y");
		next if ($pos[0] ne $main::chrom && $main::chrom ne "1-24");
		$Weight_prob{"$pos[0]:$pos[1]"}=$row[1];
	}

	close (PolyPh);


	# read functional effects

	my %Func_eff;

	open (EFFECTS, "$INSTALL_PATH/VTtest/functional_effects.txt");

	while (<EFFECTS>){
		my @row=split(/\s+/,$_);
		$row[0] =~ s/chr//;
		$row[0]=23 if ($row[0] eq "X");
		$row[0]=24 if ($row[0] eq "Y");
		next if ($row[0] ne $main::chrom && $main::chrom ne "1-24");
		$Func_eff{"$row[0]:$row[1]"}=$row[2];
	}

	close (EFFECTS);


	# read and hash VCF genotypes
	
	my $missing_function = 0;

	open (UNKNOWN, ">>unknown-functions.txt");
	
	my @Sample_IDs;
	
	my %alleles = ();
	my @hash_minor = ();
	
	open (INPUT,"$input") or die "can't open input: $input\n";

	while(<INPUT>){
	
	    chomp;

	   	if ($_ =~ /^\#/){
		   	if ($_ =~ /^\#\#/){
				next;
			} else {
				my @vcf_header = split (/\s+/,$_);
				@Sample_IDs = @vcf_header[9..$#vcf_header];
				next;
		   	}
		}


		my @array = split(/\s+/,$_);


		next if ($array[4] =~ /,/); 
		# currently doesn't handle more than two alleles
    
		$array[0]=23 if ($array[0] eq "X");
		$array[0]=24 if ($array[0] eq "Y");
	
		next if ($array[0] ne $main::chrom && $main::chrom ne "1-24");
	
		# Filter by single point fisher's exact p-value
		next if (exists $main::fisher_hash{"$array[0]:$array[1]"});

		# filter by allele frequency:
		#next if ($main::Minor_allele{"array[0]:$array[1]"} >= 0.05 
		#     || $main::Minor_allele{"array[0]:$array[1]"} <= 0.95);
		
		# Extract functional consequence and polyphen probs
		my $VT_Weight;
		
		if ($array[7] =~ /CSQ/){
			
			$Func_eff{"$array[0]:$array[1]"}=$array[7];
			
			if ($array[7] =~ /NON_SYNONYMOUS_CODING/ && $array[7] =~ /PolyPhen/){
				
				my @row=split(/:/,$array[7]);
			
				foreach my $item (@row){
									
					if ($item =~ /PolyPhen/){
						
						#print "$item\n";
						$item =~ tr/0-9\.0-9//cd;
						my $weight = `echo \"$item\" | $INSTALL_PATH/VTtest/pph2vt.pl`;
						#print "$item\t$weight\n";
						chomp($weight);
						if (defined $VT_Weight){
							if ($VT_Weight<$weight){
								$VT_Weight=$weight;
							}
						} else {
							$VT_Weight=$weight;
						};
					
					}
					
				}
				
			}
			
		}
		
		
		if (!defined $VT_Weight && exists $Weight_prob{"$array[0]:$array[1]"}){		
			$VT_Weight = $Weight_prob{"$array[0]:$array[1]"};
		}
	

		#
		
		my $missing_function = 0;

		if ($main::Minor_freq{$SNP} eq "." || $main::Minor_freq{$SNP} >= 0.01){ # common snps		

			if (defined $Func_eff{$SNP}) {
						
				if ($Func_eff{$SNP} =~/STOP_GAINED/ || $Func_eff{$SNP} =~/FRAMESHIFT/ || $Func_eff{$SNP} =~/SPLICE_SITE/) {print OUTPUT "$array[0]\t$array[1]\t0.973985";} else {print OUTPUT "$array[0]\t$array[1]\t0.373657";};
						
			} else { 
						
				print OUTPUT "$array[0]\t$array[1]\t0.373657";
							
		 		print UNKNOWN "$array[0]\t$array[1]\t$array[1]\t$array[3]/$array[4]\t+\n";
				$missing_function = 1;
						
			}				
					
		} elsif (defined $VT_Weight){ # rare snps
				
			print OUTPUT "$array[0]\t$array[1]\t".$VT_Weight;
				
		} else {
						
			print OUTPUT "$array[0]\t$array[1]\t0.373657";
						
		}


		# Identify the order of INFO fields
		my @info_fields = split (/:/, $array[8]);
		
		my ( $GT_index ) = grep { $info_fields[$_] eq "GT" } 0..$#info_fields;
		my ( $GQ_index ) = grep { $info_fields[$_] eq "GQ" } 0..$#info_fields;
		my ( $DP_index ) = grep { $info_fields[$_] eq "DP" } 0..$#info_fields;
		
		# determine count of minor alleles
		for (my $i = 9; $i <= $#array; $i++){

			my $minor_count = 0;
	
			my @info = split(/:/,$array[$i]);
	 
 		    my $first_allele = substr $info[$GT_index], 0, 1;
	 	 	my $second_allele = substr $info[$GT_index], -1, 1;
	 
		 	if ($main::Minor_allele{"$array[0]:$array[1]"} eq $array[4]){

				$minor_count++ if ($first_allele eq "1");
				$minor_count++ if ($second_allele eq "1");
			  	#$minor_count = 0 if (defined $info[2] && $info[2] < 8); # filter by depth 
		 	
		 	} elsif ($main::Minor_allele{"$array[0]:$array[1]"} eq $array[3]){
			
				$minor_count++ if ($first_allele eq "0");
				$minor_count++ if ($second_allele eq "0");
			  	#$minor_count = 0 if (defined $info[2] && $info[2] < 8); # filter by depth 
		 	
		 	} #else {
		 	
		 		#die "Minor allele mismatch for SNP: chr$array[0] $array[1]..\n";
		 	
		 	#}
	 	
		 	print OUTPUT "\t$minor_count";

		}
		
		print OUTPUT "\n";

	}

	close (INPUT);


	if ($missing_function == 1){
		print "Functional effect of a number of SNPs is unknown:  ./unknown-functions.txt\n";
		print "Retrieve their function using SNP Effect Predictor and re-run the pipeline..\n";
		print "See: http://www.ensembl.org/Homo_sapiens/UserData/UploadVariations\n";
		
		print "\nPress \"I\" to ignore and use pph2_prob=0.5 for them...\n";
		my $key = "I"; #<>;	
		chomp($key); 	die if ($key !~ m/I/i);
	}
	
	close (UNKNOWN);
	
	close(OUTPUT);


	## create individual weight and geno files

	foreach my $gene (keys %main::region_start){
    
    	my $chr = $main::region_chr{$gene};
    	
    	$chr = 23 if ($chr eq "X");
        $chr = 24 if ($chr eq "Y");


		my $geno_file = "$main::VT_geno_dir/".$gene.".cgeno.dat";
		open (GENO_OUTPUT,">$geno_file") or die "can't open output: $geno_file\n";
		
		my $weight_file = "$main::VT_weight_dir/".$gene.".polyph.dat";
		open (WEIGHT_OUTPUT,">$weight_file") or die "can't open output: $weight_file\n";
	
		my @variants = `awk '\$1==$chr && \$2>=$main::region_start{$gene} && \$2<=$main::region_end{$gene}' $vt_format_file`;

		foreach my $var (@variants){
			
			my @row=split(/\s+/,$var);
			
			print WEIGHT_OUTPUT "$row[1] $row[2]\n";
			
			for my $i (0 .. $#Sample_IDs){
			
				if (exists $main::Pheno_dicho{"$Sample_IDs[$i]"} 
						|| exists $main::Pheno_quant{"$Sample_IDs[$i]"}){
							
					#if ($row[$i+3] ne 0){
						print GENO_OUTPUT $Sample_IDs[$i]." $row[1] $row[$i+3]\n";
					#};
							
				}
	
			}	
			
		}
		
	    close (WEIGHT_OUTPUT);		    
	    close (GENO_OUTPUT);


	}

}

1;

