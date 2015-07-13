#!/usr/bin/perl

#########################################################################################################
# This script compares SNPs from Guzerá detected in this study with SNPs from Gir detected by Lee et al., 2013. 
# You may distribute this script under the same terms as perl itself
# Francislon Silva <francislon at cpqrr dot fiocruz dot br>
# GGBC Fiocruz-MG - 22_10_2014
#########################################################################################################

use strict;
use Getopt::Long;

#########################################################################################################
# Declaring variables
#########################################################################################################

my $file_gir;           # file with snps and indels from gir
my $file_snp_guzera;    # file snps from guzera
my $file_indel_guzera;     # file indels from guzera
my $file_output_shared_snp;          # output file name to export all snps shared between gir and guzera
my $file_output_shared_indels;          # output file name to export all indels shared between gir and guzera
my $file_output_excl_snp_gir;          # output file name to export all snps exclusives from gir
my $file_output_excl_indel_gir;          # output file name to export all indels exclusives from gir

my $file_output_excl_snp_guz;          # output file name to export all snps exclusives from guzera
my $file_output_excl_indel_guz;          # output file name to export all indels exclusives from gir

my $count_lines = 0;      # store the quantity of lines from each file
my %hash_count  = ();     # hash to manage the reading progress


my %hash_indel_gir=(); # hash to store all indels informations from gir
my %hash_snp_gir=(); # hash to store all snps informations from gir

my $fh_snps_shared; # filehandle to export all snps shared between gir and guzera
my $fh_snps_guzera;  # filehandle to export all snps exclusives from guzera
my $fh_snps_gir;  # file handle to export all snps exclusives from gir

my $fh_indels_shared; # filehandle to export all indels shared between gir and guzera
my $fh_indels_guzera; # filehandle to export all indels from guzera
my $fh_indels_gir; # filehandle to export all indels from gir


my $total_snp_shared=0; # total of snps shared between gir and guzera
my $total_indels_shared=0; # total of indels shared between gir and guzera

my $total_indels_gir=0; # total of indels from gir
my $total_snp_gir=0; # total of snps from gir

my $total_indels_guz=0; # total of indels from guzera
my $total_snp_guz=0; # total of snps from guzera

my $header_gir = ""; # header for output

my $help; # flag to know if help is requested

#########################################################################################################
# Receiving parameters
#########################################################################################################
GetOptions(
			'g=s'    => \$file_gir,
			's=s'    => \$file_snp_guzera,
			'z=s'    => \$file_indel_guzera,
			'oss=s'    => \$file_output_shared_snp,
			'osi=s'    => \$file_output_shared_indels,
			'oeig=s'    => \$file_output_excl_indel_gir,
			'oesg=s'    => \$file_output_excl_snp_gir,
			'oeiz=s'    => \$file_output_excl_indel_guz,
			'oesz=s'    => \$file_output_excl_snp_guz,
			'help|?' => \$help
);

$| = 1; # forces a flush right away and after every write or print on the currently selected output channel

#########################################################################################################
# Calling the main function
#########################################################################################################
main();


########################################################################
# Main function of the program
########################################################################
sub main {

	validate_parameters(); # Calling the function to validate the parameters
	
	print "Reading VCF file for Gir...\n";
	read_file_gir();    # reading bioscope input file and generating output
	
	open_files_for_write(); # opening all files to write
	
	print "\nReading tab file for Guzera...\n";
	read_snp_guzera(); # Reading guzera's snps from file
	
	print "\nReading VCF file for Guzera...\n";
	read_indel_guzera(); # reading VCF file with indels from guzera
	close_files_for_write(); # closing all files
	print "\n-------------------------------------------------\n";
	print "Total SNPs shared: $total_snp_shared\n";
	print "Total INDELs shared: $total_indels_shared\n";
	
	print "Gir exclusives: \n";
	print "\tSNPs: $total_snp_gir\n";
	print "\tINDELs: $total_indels_gir\n";
	
	print "Guzera exclusives: \n";
	print "\tSNPs: $total_snp_guz\n";
	print "\tINDELs: $total_indels_guz\n";

	print "\nProcess finished!\n";

}

#########################################################################################################
# This function opens all files to write
#########################################################################################################
sub open_files_for_write{
	open($fh_snps_shared, ">".$file_output_shared_snp);
	open($fh_indels_shared, ">".$file_output_shared_indels);
	
	open($fh_indels_gir, ">".$file_output_excl_indel_gir);
	open($fh_snps_gir, ">".$file_output_excl_snp_gir);
	
	open($fh_indels_guzera, ">".$file_output_excl_indel_guz);
	open($fh_snps_guzera, ">".$file_output_excl_snp_guz);
}

#########################################################################################################
# This function closes all files to write
#########################################################################################################
sub close_files_for_write{
	close($fh_snps_shared);
	close($fh_indels_shared);
	
	close($fh_indels_gir);
	close($fh_snps_gir);
	
	close($fh_indels_guzera);
	close($fh_snps_guzera);
}

#########################################################################################################
# This function reads a VCF file containing indels from Guzera 
#########################################################################################################
sub read_indel_guzera{
	
	open( IN, $file_indel_guzera );
	
	my $header = "";
	my $print_header = 0;
	
	my $count                     = 0; # count the lines read
	my $wrap_line                 = 0; # flag to determine when break a line	
	%hash_count = (); # hash to manage the percentage of reading progress
	$count_lines = `wc -l $file_indel_guzera`; # necessary to show the reading progress
	
	print "Progress: ";
	while (<IN>) {
		chomp;
		unless (/^(\s|\t)*$/) { # eliminating blank lines
			unless (/^#/) { # reading 
				if(!$print_header){
					print $fh_indels_shared $header;
					print $fh_indels_guzera $header;
					$print_header = 1;
				}
				my $line = $_;
				#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Unknown
				#  0	1	2	3	4	5		6		7		8		9
				my @columns_inner = split /\t/;    # spliting line by tab

				my $info = $columns_inner[7];
				my $chr  = $columns_inner[0];
				$chr =~ s/Chr(\d+)/$1/g;
				my $pos  = $columns_inner[1];

				if ( $info =~ /^INDEL;/ ) {
					if(exists $hash_indel_gir{$chr}{$pos}){
						$total_indels_shared++;
						print $fh_indels_shared $line."\n";
						delete $hash_indel_gir{$chr}{$pos};
					}else{
						$total_indels_guz++;
						print $fh_indels_guzera $line."\n";
					}
				}
				
				# the following lines is just to calculate and printing the percentage of progress
				my $percent_lines = int( ( $count++ * 100 ) / $count_lines );
				unless ( exists $hash_count{$percent_lines} ) {
					print " | $percent_lines%";
					$wrap_line++;
					$hash_count{$percent_lines} = $percent_lines;
					if ( $percent_lines > 0 and $wrap_line % 10 == 0 ) {
						print "\n";
					}
				}
				
							
			}
			else {
				$header .= $_;
				$header .= "\n";
				$count_lines--;
			}
		}
		else {
			$count_lines--;
		}
	}
	print "\n";
	close(IN);
	
	# exporting gir exclusives
	print "\nExporting INDELs exclusives for Gir...\n";
	print $fh_indels_gir $header_gir;
	
	my @chromosomes = sort {$a<=>$b} keys %hash_indel_gir;
	for(@chromosomes){
		my $chr = $_;
		my @positions = sort {$a<=>$b} keys %{$hash_indel_gir{$_}};
		$total_indels_gir += scalar @positions;
		for(@positions){
			print $fh_indels_gir $hash_indel_gir{$chr}{$_}."\n";	
		}
	}
	
}

#########################################################################################################
# This function reads a TAB file containing snps from guzera
#########################################################################################################
sub read_snp_guzera{
	
	
	# opening the input file and the output file
	open( IN, $file_snp_guzera );

	# reading the header
	my $header = <IN>;
	my @lines_filtered_inner = ();
	my $count                     = 0;
	my $wrap_line                 = 0;
	
	%hash_count = ();
	$count_lines = `wc -l $file_snp_guzera`;


	print $fh_snps_shared $header;
	print $fh_snps_guzera $header;

	print "Progress: ";
	while (<IN>) {
		chomp;
		unless (/^(\s|\t)*$/) {
			my $line = $_;
# Seqid	Source	Type	Start	End	Score	Strand	Phase	genotype	reference	coverage	refAlleleCounts	refAlleleStarts	refAlleleMeanQV	novelAlleleCounts	novelAlleleStarts	novelAlleleMeanQV	mostAbundantAlleleDiColor2	secondAbundantAlleleDiColor3	het	geneID	exonID	rsID	functionCode
#	0	1		2		3		4	5		6		7		8			9			10			11				12				13				14					15					16					17							18								19	20		21		22		23
			my @columns_inner = split /\t|\s+/;    # spliting line by tab
			my $chr = $columns_inner[0];
			$chr =~ s/^Chr(\d+)/$1/gi;
			my $pos = $columns_inner[3];
			
			if(exists $hash_snp_gir{$chr}{$pos}){
				print $fh_snps_shared $line."\n";
				delete $hash_snp_gir{$chr}{$pos};
				$total_snp_shared++;
			}else{
				print $fh_snps_guzera $line."\n";
				$total_snp_guz++;
			}
			
			# the following lines is just to calculate and printing the percentage of progress
				my $percent_lines = int( ( $count++ * 100 ) / $count_lines );
				unless ( exists $hash_count{$percent_lines} ) {
					print " | $percent_lines%";
					$wrap_line++;
					$hash_count{$percent_lines} = $percent_lines;
					if ( $percent_lines > 0 and $wrap_line % 10 == 0 ) {
						print "\n";
					}
				}
			
		}
		else {
			$count_lines--;
		}
	}
	print "\n";
	close(IN);
	
	print "\nExporting SNPs exclusives for Gir...\n";
	# exporting gir exclusives
	print $fh_snps_gir $header_gir;
	my @chromosomes = sort {$a<=>$b} keys %hash_snp_gir;
	for(@chromosomes){
		my $chr = $_;
		my @positions = sort {$a<=>$b} keys %{$hash_snp_gir{$_}};
		$total_snp_gir += scalar @positions;
		for(@positions){
			print $fh_snps_gir $hash_snp_gir{$chr}{$_}."\n";	
		}
	}
	
}

#########################################################################################################
# This function reads a VCF file containing snps and indels from gir
#########################################################################################################
sub read_file_gir {
	open( IN, $file_gir );

	my $count                     = 0;
	my $wrap_line                 = 0;
	
	%hash_count = ();
	$count_lines = `wc -l $file_gir`;

	print "Progress: ";
	while (<IN>) {
		chomp;
		unless (/^(\s|\t)*$/) {
			unless (/^#/) {
				my $line = $_;
				#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Unknown
				#  0	1	2	3	4	5		6		7		8		9
				my @columns_inner = split /\t/;    # spliting line by tab

				my $info = $columns_inner[7];
				my $chr  = $columns_inner[0];
				my $pos  = $columns_inner[1];
				
				if($info =~ /VRT=1/){ # SNP
					$hash_snp_gir{$chr}{$pos} = $line;
				}else{ # INDEL
					$hash_indel_gir{$chr}{$pos} = $line;
				}
				# the following lines is just to calculate and printing the percentage of progress
				my $percent_lines = int( ( $count++ * 100 ) / $count_lines );
				unless ( exists $hash_count{$percent_lines} ) {
					print " | $percent_lines%";
					$wrap_line++;
					$hash_count{$percent_lines} = $percent_lines;
					if ( $percent_lines > 0 and $wrap_line % 10 == 0 ) {
						print "\n";
					}
				}
				
				
			}
			else {
				$header_gir .= $_;
				$header_gir .= "\n";
				$count_lines--;
			}
		}
		else {
			$count_lines--;
		}
	}
	print "\n";
	close(IN);
}

#########################################################################################################
# This function validates the parameters provided from user to the script
#########################################################################################################
sub validate_parameters {
	my $allExists  = 1;
	my $fileExists = 1;

	if ( defined $help ) {
		print usage();
		exit 0;
	}

	unless ( defined $file_gir ) {
		$allExists = 0;
	}

	unless ( defined $file_snp_guzera ) {
		$allExists = 0;
	}

	unless ( defined $file_indel_guzera ) {
		$allExists = 0;
	}
	
	unless ( defined $file_output_excl_indel_gir ) {
		$allExists = 0;
	}
	
	unless ( defined $file_output_excl_indel_guz ) {
		$allExists = 0;
	}
	
	unless ( defined $file_output_excl_snp_gir ) {
		$allExists = 0;
	}
	
	unless ( defined $file_output_excl_snp_guz ) {
		$allExists = 0;
	}
	
	unless ( defined $file_output_shared_indels ) {
		$allExists = 0;
	}
	
	unless ( defined $file_output_shared_snp ) {
		$allExists = 0;
	}
	
	

	if ($allExists) {
		unless ( -e $file_gir ) {
			print STDERR "$file_gir doesn't exists.\n";
			$fileExists = 0;
		}
		unless ( -e $file_snp_guzera ) {
			print STDERR "$file_snp_guzera doesn't exists.\n";
			$fileExists = 0;
		}
		unless ( -e $file_indel_guzera ) {
			print STDERR "$file_indel_guzera doesn't exists.\n";
			$fileExists = 0;
		}
	}

	unless ($allExists) {
		print usage();
		exit 0;
	}

	unless ($fileExists) {
		print STDERR "Program execution aborted.\n";
		exit 0;
	}

}

#########################################################################################################
# This function return the text explaining how to use the script
#########################################################################################################
sub usage {
	my $usage = <<FOO;
Usage:
	perl $0 -g gir_indels_and_snp -z guzera_indels -s guzera_snps -oss shared_snps -osi shared_indels -oeig exclusives_indels_gir -oesg exclusives_snps_gir -oeiz exclusives_indels_guzera -oesz exclusives_snps_guzera 
	
	gir_indels_and_snp = VCF file containing snps and indels from gir
	guzera_indels = VCF file containing indels from guzera
	guzera_snps = TAB file containing snps from guzera
		
	The parameters started with '-o' are output files.	
		
FOO
	return $usage;
}
