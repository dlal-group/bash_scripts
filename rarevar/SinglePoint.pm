package SinglePoint;

sub CreateFisherTable {

	my $input = shift; #input file in VCF format
	
	#print "inside SinglePoint: $input\n";

	#my @path = split(/\//,$input);
	my $output = $input; #pop (@path);
	$output =~ s/\.recode\.vcf/\.fisher/;


	print "Creating contingency tables for genotypes...\n\n";

	my @Sample_IDs;

	open (OUTPUT, ">$output") or die "can't open input: $output\n";
	open (INPUT,"$input") or die "can't open input: $input\n";

	print OUTPUT "chr\tpos\tcontrols_ref\tcontrols_alt\tcases_ref\tcases_alt\n";

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

		my @array = split(/\t/,$_);

		$array[0] = 23 if ($array[0] eq "X");
		$array[0] = 24 if ($array[0] eq "Y");
		
		next if ($array[0] ne $main::chrom && $main::chrom ne "1-24");
		
		
		# Skip loci with more than two alleles
		next if ($array[4] =~ /,/); 


		# Identify the order of INFO fields
		my @info_fields = split (/:/, $array[8]);
		
		my ( $GT_index ) = grep { $info_fields[$_] eq "GT" } 0..$#info_fields;
		my ( $GQ_index ) = grep { $info_fields[$_] eq "GQ" } 0..$#info_fields;
		my ( $DP_index ) = grep { $info_fields[$_] eq "DP" } 0..$#info_fields;
		
		#print "$GT_index\t$GQ_index\t$DP_index\n";

		my $Cases_alt_count = 0;
		my $Cases_ref_count = 0;
	
		my $Controls_alt_count = 0;
		my $Controls_ref_count = 0;

	
		for (my $i = 9; $i <= $#array; $i++){
	
			my $my_id = $Sample_IDs[$i-9];
		
			next if (!exists $main::Pheno_dicho{$my_id});
		
			#print "$i $my_id ".$main::Pheno_dicho{$my_id}."\n";
			
			my @info = split(/:/,$array[$i]);
	 
			if (defined $info[$GT_index] && $info[$GT_index] ne "./."){

				my $first_allele = substr $info[$GT_index], 0, 1;
 	 			my $second_allele = substr $info[$GT_index], -1, 1;	
			
				if (exists $main::Pheno_controls{$my_id}) {	
					#print "Cases: $Cases_ref_count $Cases_alt_count\n";				

					$Cases_ref_count++ if ($first_allele eq "0");
					$Cases_alt_count++ if ($first_allele eq "1");

					$Cases_ref_count++ if ($second_allele eq "0");
					$Cases_alt_count++ if ($second_allele eq "1");

				} elsif (exists $main::Pheno_cases{$my_id}) {		
					#print "Controls: $Controls_ref_count $Controls_alt_count\n";
				
					$Controls_ref_count++ if ($first_allele eq "0");
					$Controls_alt_count++ if ($first_allele eq "1");

					$Controls_ref_count++ if ($second_allele eq "0");
					$Controls_alt_count++ if ($second_allele eq "1");				
				
				}
			
			}

		}
	
		print OUTPUT 
		"$array[0]\t$array[1]\t$Cases_ref_count\t$Cases_alt_count\t$Controls_ref_count\t$Controls_alt_count\n";
	
	}

	close (INPUT);
	close (OUTPUT);

	return ($output);
}

1;
