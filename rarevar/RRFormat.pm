package RRFormat;

sub FormatChange{


	my $input = shift;    #input VCF format

	#my @path = split(/\//,$input);
	my $output = $input; #pop (@path);
	$output =~ s/\.recode\.vcf/\.RR_format\.txt/;


	open(INPUT,"$input") or die "can't open input: $input\n";
	open(OUTPUT,">$output") or die "can't open output: $output\n";

	$num_lines = 0;

	while(<INPUT>){

	    chomp;
	    #if ($num_lines > 100){last;}
	    if ($_ =~ /^#/){next;}
	    else
	      {
		@array = split(/\t/,$_);
		#if ($array[6] != 0){next;}
		$last_index = $#array;
		if ($array[5] eq '.' || $array[5]>100){$qual = 100;}#if not qual value set to 100
		else                 {$qual = $array[5];}

		$array[0] = 23 if ($array[0] eq "X");
		$array[0] = 24 if ($array[0] eq "Y");
		
		next if ($array[0] ne $main::chrom && $main::chrom ne "1-24");

		# Filter by single point fisher's exact p-value
		next if (exists $main::fisher_hash{"$array[0]:$array[1]"});

		print OUTPUT "$array[0]\t$array[1]\t$array[2]\t$array[3]\t$array[4]";



		$geno_index = '';
		$qual_index = '';#if not qual index set qual to 100 everywhere
		@format_array = split(/:/,$array[8]);
		$format_last_index = $#format_array;
		for ($i = 0; $i <= $format_last_index; $i++)
		  {
		    if ($format_array[$i] eq 'GT'){$geno_index = $i;}
		  }
		if ($geno_index !~ /^\d+$/){die "There isn't a genotype field to convert\n";}

		for ($i = 9; $i <= $last_index; $i++)
		  {
		    @genotype_info_array = split(/:/,$array[$i]);
		    $genotype            = $genotype_info_array[$geno_index];

			
			print OUTPUT "\t";

		    if ($genotype =~ /\|/)
		      {
			($first_allele,$second_allele) = split(/\|/,$genotype);
			print OUTPUT $first_allele+$second_allele;
		      }
		    elsif ($genotype =~ /\//)
		      {
			($first_allele,$second_allele) = split(/\//,$genotype);
			print OUTPUT $first_allele+$second_allele;	
		      }
		    elsif ($genotype =~ /\\/)
		      {
			($first_allele,$second_allele) = split(/\\/,$genotype);
			print OUTPUT $first_allele+$second_allele;
		      }
		    else
		      {
			$first_allele = ".";
			$second_allele = ".";
			print OUTPUT ".";
		      }  
		  }#end foreach individual
		print OUTPUT "\n";
		$num_lines += 1;
	#	if ($num_lines == 1){$first_bp_pos = $array[1];}
	#	else
	#	  {
	#	    $diff = $array[1] - $first_bp_pos;
	#	    if ($diff > 150000){last;}
	#	  }
	      }#end else data line
	  }#end while input
	close(INPUT);
	close(OUTPUT);

	return ($output);

}

1;
