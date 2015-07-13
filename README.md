# Guzera repository

This repository aims to store all source code used in Guzera Project. For now we have one script available:

* [findGenesInQTLRegions.pl](/findGenesInQTLRegions.pl/)
* [compare_snps_indels.pl](/compare_snps_indels.pl/)

## findGenesInQTLRegions.pl

### Prerequisite

1. Perl interpreter installed on your operating system.

### Running the script

To run the script just type on terminal:

$ perl findGenesInQTLRegions.pl -g gff_file -q qtl_file -d dbsnp_file -o output_file [-h genes_list]

## compare_snps_indels.pl

### Prerequisite

1. Perl interpreter installed on your operating system.

### Running the script

To run the script just type on terminal:

$ perl compare_snps_indels.pl -g gir_indels_and_snp -z guzera_indels -s guzera_snps -oss shared_snps -osi shared_indels -oeig exclusives_indels_gir -oesg exclusives_snps_gir -oeiz exclusives_indels_guzera -oesz exclusives_snps_guzera 
	
	gir_indels_and_snp = VCF file containing snps and indels from gir
	guzera_indels = VCF file containing indels from guzera
	guzera_snps = TAB file containing snps from guzera


This documentation is under construction. Please send an email to francislon@gmail.com if you have any question.
