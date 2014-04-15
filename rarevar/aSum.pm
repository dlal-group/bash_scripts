package aSum;

sub vcf2geno {

	my $input = shift; #input file in VCF format


	# read and hash VCF genotypes

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


		# Identify the order of INFO fields
		my @info_fields = split (/:/, $array[8]);
		
		my ( $GT_index ) = grep { $info_fields[$_] eq "GT" } 0..$#info_fields;
		my ( $GQ_index ) = grep { $info_fields[$_] eq "GQ" } 0..$#info_fields;
		my ( $DP_index ) = grep { $info_fields[$_] eq "DP" } 0..$#info_fields;
		
		my $total_count = 0;
		
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
	 	
		 	$hash_minor[$i-9]{"$array[0]:$array[1]"} = $minor_count;
		 	
		 	$total_count += $minor_count;

		}

		$alleles{"$array[0]:$array[1]"}="$array[3]/$array[4]" if ($total_count!=0); #only polymorphic SNPs
		
	}

	close (INPUT);



	### print out genotype matrices

	my @all_snps = keys %alleles;
	

	for my $chr_id (1 .. 24){

		next if ($chr_id ne $main::chrom && $main::chrom ne "1-24");
	
		my $pat = "^$chr_id:";
	
	    my @chr_snps = grep(/$pat/, @all_snps);
       
        #print "@chr_snps";
       
	    foreach my $gene (keys %main::region_start){
    
                $chr_id = "X" if ($chr_id eq "23");
                $chr_id = "Y" if ($chr_id eq "24");
	
	    	next if ($main::region_chr{$gene} ne $chr_id);
    		
			my $output = "$main::aSum_geno_dir/".$gene.".geno.dat";
			#print "$output..\n";
			open (GENO_OUTPUT,">$output");# or die "can't open output: $output\n"
	
	
			foreach my $SNP (@chr_snps){
				
				my @pos = split (/:/,$SNP);	
		
				if ($pos[1] >= $main::region_start{$gene} && $pos[1] <= $main::region_end{$gene}) {
						
					for my $i (0 .. $#Sample_IDs){
					
						# print minor allele counts in each regions
						if (exists $main::Pheno_dicho{"$Sample_IDs[i]"})
							{print GENO_OUTPUT $hash_minor[$i]{$SNP}." ";}
					
					}
					
					print GENO_OUTPUT "\n";
	    
		    	}
			
			}
			
			close (GENO_OUTPUT);
	
		}
	
	}
	

}

1;
