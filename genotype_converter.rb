#!/usr/bin/ruby
##Script to convert genotypes in the correct format using a reference file
#we can use also for genotype discordance calculation

#open the ref file
#FIXME:set a variable name for ref file...
ref_file=File.open("ref_discordance_table.txt")

ref_file.each do |line|
	puts line.split("\t")[1]
end
