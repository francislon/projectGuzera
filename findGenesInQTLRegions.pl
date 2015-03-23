#!/usr/bin/perl

#########################################################################################################
# Recover the genes associated with traits of interest using the flank markers information available on Cattle QTL in the Animal QTLdb.
# You may distribute this script under the same terms as perl itself
# Francislon Silva <francislon at cpqrr dot fiocruz dot br> and Izinara Rosse <izinara dot rossi at gmail dot com> 
# GGBC Fiocruz-MG - 22_10_2014
#########################################################################################################

use strict;
use Getopt::Long;

#########################################################################################################
# Declaring variables
#########################################################################################################
use constant{
	NAME => 'name',
	POS_START => 'start',
	POS_END => 'END',
	QTL=>'qtl',
	CATEGORY_QTL=>'cat_qtl',
	INFO	=>'info',
	DBSNP_ID	=>'dbsnpid',
	FILE	=> 'file'
};

my $file_genes_list; # export genes that exists in QTL regions
my $file_dbsnp; # dbsnp database file
my $file_gff; # all genes annotated
my $file_qtl; # qtl database file
my $file_output; # output file with SNPs found in any gene described in a GFF file
my $file_output_no_rs; # output file with 'FlankersMarkers' (from QTL database file) not found in a dbsnp file
my $file_output_gene_qtl; # output file with genes from any output exported in this script
my $file_output_leftovers; # output file with SNPs around 100 kb from any gene described in a GFF file

my $count_lines_gene_list = 0; # quantity of lines in a $file_genes_list
my $count_lines_dbsnp = 0; # quantity of lines in a $file_dbsnp
my $count_lines_gff = 0; # quantity of lines in a $file_gff
my $count_lines_qtl = 0; # quantity of lines in a $file_qtl
my %hash_count = (); # hash to manage the file reading progress

my %hash_qtl=(); # hash to store all QTL information from QTL file
my %hash_dbsnp=(); # hash to store all dbSNP information from dbSNP file
my %hash_snp_position_gene=(); # hash to store snp position and gene related
my %hash_leftovers=(); # hash to store position of SNPs around 100 kb from any gene described in a GFF file

my %hash_leftovers_genes=(); # hash to store information about genes with SNPs 100 kb around the gene

my %hash_leftovers_inserted=();# hash to control what genes can be exported in the output

my %hash_genes_list=(); # hash to store genes from $file_genes_list [only to facilitate the search]
my @array_genes_list=(); # array to store genes from $file_genes_list [to keep the order of the file]
my %hash_genes_in_files=(); # hash to store genes that will be exported if the $file_genes_list is provided

my $help; # flag to know if help is requested

#########################################################################################################
# Receiving parameters
#########################################################################################################
GetOptions('d=s' => \$file_dbsnp,
			'g=s' => \$file_gff,
			'q=s' => \$file_qtl,
			'h=s' => \$file_genes_list,
			'o=s' => \$file_output,
			'help|?' => \$help
);

$| = 1; # forces a flush right away and after every write or print on the currently selected output channel

#########################################################################################################
# Calling the main function
#########################################################################################################
main();

#########################################################################################################
# The main function
#########################################################################################################
sub main{
	
	validate_parameters(); # This function validates if the parameters are correct
	
	pre_processing(); # This function counts the quantity of lines from each file
	
	# if $file_genes_list is provided then we will read the file
	if(defined $file_genes_list){
		print "Reading $file_genes_list...\n";
		$file_output_gene_qtl = $file_output.".gene_qtl";
		read_genes_list();
	}
	
	print "Reading $file_qtl...\n";
	read_qtl();
	print "Reading $file_dbsnp...\n";
	read_dbsnp();
	
	print "Reading $file_gff...\n";
	read_gff();
	
	if(defined $file_genes_list){
		export_gene_qtl();
	}
	
	print "\n\nFiles [$file_output], [$file_output".".flank_no_rs]";
	if(defined $file_genes_list){
		print ", [$file_output_gene_qtl] and [$file_output".".remaining] generated...\n\n"
	}
	print "Process finished!\n";
	
}

#########################################################################################################
# This function counts the quantity of lines from each file
#########################################################################################################
sub pre_processing{
	if(defined $file_genes_list){
		$count_lines_gene_list = `wc -l $file_genes_list`;	
	}
	$count_lines_dbsnp = `wc -l $file_dbsnp`;
	$count_lines_gff = `wc -l $file_gff`;
	$count_lines_qtl = `wc -l $file_qtl`;
}

#########################################################################################################
# This function reads a file with one gene per line
# If the gene list file is provided we will export all genes in that file who was exported 
# in another exported file
#########################################################################################################
sub read_genes_list{
	
	my $count = 0;
	my $wrap_line = 0;
	%hash_count = ();
	
	open(IN, $file_genes_list);
	
	while(<IN>){
		chomp;
		my $line = $_;
		unless(/^(\s|\t)*$/){
			my @columns = split /\s+/;
			
			$columns[0] =~ s/^>//g;
			$columns[0] =~ s/\s//g;
			$hash_genes_list{$columns[0]} = 1;
			push(@array_genes_list, $columns[0]);
			
		} #end unless(/^(\s|\t)*$/
		
		# the following lines are just to calculate and printing the percentage of progress
		my $percent_lines = int(($count++ * 100)/$count_lines_gene_list);		
		unless(exists $hash_count{$percent_lines}){
			print " | $percent_lines%";
			$wrap_line++;
			$hash_count{$percent_lines} = $percent_lines;
			if( ($percent_lines > 0 and $wrap_line % 10 == 0) or $percent_lines==100 ){
				print "\n";
			}
		}
	}
	
	unless(exists $hash_count{100}){
		print " | 100%\n";
	}
	close(IN);
}


#########################################################################################################
# This function reads a QTL file.
# We are looking for genes in QTL regions that we can verify if there are SNPs in those genes
#########################################################################################################
sub read_qtl{
	
	my $count = 0;
	my $wrap_line = 0;
	%hash_count = ();
	
	open(IN, $file_qtl);
	$file_output_no_rs = $file_output.".flank_no_rs";
	# the file $file_output_no_rs will export all QTL regions which have at least one gene starting with no 'rs'
	open(OUT, ">".$file_output_no_rs);
	
	while(<IN>){
		chomp;
		my $line = $_;
		unless(/^(\s|\t)*$/){
			unless(/^#/){
				my @columns = split /\t/;
				my $chr = $columns[0];
				$chr =~ s/^Chr\.(\w+)$/$1/g;
				my $qtl = $columns[1];
				my $category_qtl = $columns[2];
				my $start = $columns[3];
				my $end = $columns[4];
				my $info_column = $columns[8];
				
				my %hash = (POS_START=>$start, POS_END=>$end,
				QTL=>$qtl,INFO=>$info_column);
				
				$chr = lc($chr);
				if($chr == 'x'){ # chromosome X will be called 30 to facilitate the comparison between files
					$chr = 30;
				}
				
				
				my $flankers_markers;
				if($info_column =~ /^.+;FlankMarkers\=([^;]+);.+$/){ # the gene information in each QTL region is located at the property 'FlankMarkers'
					$flankers_markers = $1;
					
					my @genes_flankers = split(/,/, $flankers_markers); # each 'FlankMarkers' can have more than one gene separated by comma
					
					my $has_no_rs = 0;
					for(@genes_flankers){
						my $gene = $_;
						if(/^rs.+$/){ # we are looking for genes with rs*, which means that the gene has SNPs in dbSNP
							if(exists $hash_qtl{$chr}{$gene}{$category_qtl}){
								my $array = $hash_qtl{$chr}{$gene}{$category_qtl};
								push(@$array, \%hash);
							}else{
								my @array=(\%hash);
								$hash_qtl{$chr}{$gene}{$category_qtl} = \@array;
							}
						}else{
							if(exists $hash_genes_list{$gene}){
								my %hash = (NAME => $gene, CATEGORY_QTL=>$category_qtl, DBSNP_ID=>'not found', FILE=>$file_output_no_rs);
								if(exists $hash_genes_in_files{$gene}){
									my $array = $hash_genes_in_files{$gene};
									push(@$array, \%hash);
								}else{
									my @array=(\%hash);
									$hash_genes_in_files{$gene} = \@array;
								}
							}
							$has_no_rs = 1;
						}# end if ^rs
					}
					if($has_no_rs){ # if this QTL has at least one gene name starting with no 'rs'
						print OUT "$line\n";
					}
				}
			}
		}
		# the next lines are just to calculate and printing the percentage of progress
		my $percent_lines = int(($count++ * 100)/$count_lines_qtl);		
		unless(exists $hash_count{$percent_lines}){
			print " | $percent_lines%";
			$wrap_line++;
			$hash_count{$percent_lines} = $percent_lines;
			if( ($percent_lines > 0 and $wrap_line % 10 == 0) or $percent_lines==100 ){
				print "\n";
			}
		}
	}
	unless(exists $hash_count{100}){
		print " | 100%\n";
	}
	
	close(OUT);
	close(IN);
}

#########################################################################################################
# This function reads a GFF file.
# We will save all genes from this GFF in a hash to search if this genes have SNPs from dbSNP.
#########################################################################################################
sub read_gff{
	
	open(IN, $file_gff);
	open(OUT, ">".$file_output);
	
	print OUT sprintf("#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
	"Chr", "'Gene Name'","'Gene Start'", "'Gene End'", "'dbSNP ID'", "'dbSNP pos'", "'QTL db'", "'QTL Category'", "'QTL Start'", "'QTL End'", "'QTL Info'");
	
	my $count = 0;
	my $wrap_line = 0;
	%hash_count = ();

	my %hash_order=();
	while(my ($chr, $hash) = each(%hash_snp_position_gene)){
		my @array_pos = sort{$a<=>$b} keys %$hash;
		$hash_order{$chr}=\@array_pos;
	}
	
	my $last_end=-1;
	my $last_start=-1;
	my %hash_backup_pos=();
	my %hash_backup_genes=();
	my $c = 0;
	while(<IN>){
		chomp;
		unless(/^(\s|\t)*$/){
			
			unless(/^(#|GS)/){
				my @columns = split /\t/;
				
				$columns[2] = lc($columns[2]);
				if($columns[2] eq "gene"){
					my $last_col = $columns[-1];
					
					my @array_last_col = split(/;/, $last_col);
					
					if($array_last_col[0] =~ /^gene_id "GK(\d{6})\.\d:(\w+)".*$/){
						my $chr = int($1);
						my $gene = $2;
						if(exists $hash_order{$chr}){
							my $array_pos = $hash_order{$chr};
							my $pos = 0;
							
							if($columns[3] > $columns[4]){
								my $aux = $columns[3];
								$columns[3] = $columns[4];
								$columns[4] = $aux;
							}
							my $start = $columns[3];
							my $end = $columns[4];
							if(!exists $hash_backup_pos{$chr}){
								my @array = ();
								$hash_backup_pos{$chr} = \@array;
							}
							my $array_backup_pos = $hash_backup_pos{$chr};
							do{
								if(scalar(@$array_pos)>0){
									$pos = shift(@$array_pos);
									if($pos >= $start and $pos <= $end){
									 	print OUT sprintf("%s\t%s\t%d\t%d\t%s\t%d",
									 	 $chr, $gene,$start, $end, $hash_snp_position_gene{$chr}{$pos}, $pos);
									 	my $hash_qtl_innter = $hash_qtl{$chr}{$hash_snp_position_gene{$chr}{$pos}};
									 	while(my ($key, $value) = each(%$hash_qtl_innter)){
									 		for(@$value){
									 			print OUT sprintf("\t%s\t%s\t%d\t%d\t%s", $_->{QTL}, $key, $_->{POS_START},
									 				$_->{POS_END}, $_->{INFO});
									 		}
									 	}
									 	print OUT "\n";
									 	
									 	if(exists $hash_genes_list{$gene}){
									 		my $category_qtl = "";
									 		my $hash_qtl_innter = $hash_qtl{$chr}{$hash_snp_position_gene{$chr}{$pos}};
										 	while(my ($key, $value) = each(%$hash_qtl_innter)){
										 		if($category_qtl eq ""){
									 				$category_qtl = $key;
									 			}else{
									 				$category_qtl .= ",".$key;
									 			}									 			
										 	}
									 		
									 		
											my %hash = (NAME => $gene, CATEGORY_QTL=>$category_qtl, DBSNP_ID=>$hash_snp_position_gene{$chr}{$pos}, FILE=>$file_output);
											if(exists $hash_genes_in_files{$gene}){
												my $array = $hash_genes_in_files{$gene};
												push(@$array, \%hash);
											}else{
												my @array=(\%hash);
												$hash_genes_in_files{$gene} = \@array;
											}
										}
									 	
									 	
									 }elsif($pos <= $end){ 
									 	if(($pos >= $start-100000 and $pos <= $start) or
									 			($pos <= $end+100000 and $pos >= $end) ){
									 			if(not exists $hash_leftovers_inserted{$chr}{$pos}{$gene}){
									 				if(not exists $hash_leftovers{$chr}{$pos}){
										 				$hash_leftovers{$chr}{$pos} = 1;
										 			}
										 			
										 			my %hash_inner = (NAME=>$gene, POS_START=>$start, POS_END=>$end);
										 			if(exists $hash_leftovers_genes{$chr}{$pos}){
										 				my $array_left = $hash_leftovers_genes{$chr}{$pos};
										 				push(@$array_left, \%hash_inner);
										 			}else{
										 				my @array_left = (\%hash_inner);
										 				$hash_leftovers_genes{$chr}{$pos} = \@array_left;
										 			}
										 			push(@$array_backup_pos, $pos);	
										 			$hash_leftovers_inserted{$chr}{$pos}{$gene} = 1;
									 			}
									 		}else{
									 			push(@$array_backup_pos, $pos);
									 		}
									 }else{
									 	unshift(@$array_pos, $pos);
									 }
								}
							}while($pos <= $end and scalar(@$array_pos)>0);
							$last_end = $end;							
							$last_start = $start;
							if(not exists $hash_backup_genes{$chr}){
								my @array = ();
								$hash_backup_genes{$chr} = \@array;
							}
							my $array_genes_by_chr = $hash_backup_genes{$chr};
							my %hash_gene = (NAME=>$gene, POS_START=>$start, POS_END=>$end);
							push(@$array_genes_by_chr, \%hash_gene);
						}
					}
				}
				
			}
		}
		# the next lines are just to calculate and printing the percentage of progress
		my $percent_lines = int(($count++ * 100)/$count_lines_gff);		
		unless(exists $hash_count{$percent_lines}){
			print " | $percent_lines%";
			$wrap_line++;
			$hash_count{$percent_lines} = $percent_lines;
			if( ($percent_lines > 0 and $wrap_line % 10 == 0) or $percent_lines==100 ){
				print "\n";
			}
		}
		
	}
	unless(exists $hash_count{100}){
		print " | 100%\n";
	}
	close(IN);
	close(OUT);

	while(my ($key, $value) = each(%hash_order)){
		if(exists $hash_backup_pos{$key}){
			my $array = $hash_backup_pos{$key};
			push(@$array, @$value);
			
		}else{
			$hash_backup_pos{$key} = $value;
		}
	}
	
	while(my ($chr, $value) = each(%hash_backup_pos)){
		
		my @array_pos = sort {$a<=>$b} @$value;
		my $array_genes = $hash_backup_genes{$chr};
		
		for(@array_pos){
			my $pos = $_;
			
			my $quant_remove = 0;
			my $pos_remove=0;
			for(my $i = 0; $i < scalar @$array_genes; $i++){
				my $hash_inner = $array_genes->[$i];
				my $start = $hash_inner->{POS_START};
				my $end = $hash_inner->{POS_END};
				my $gene = $hash_inner->{NAME};
				if($pos + 100000 < $start){
					$quant_remove++;
					$pos_remove = $i;
					last;
				}elsif( ($pos >= $start-100000 and $pos <= $start) or
						($pos <= $end+100000 and $pos >= $end)){
					
					if(not exists $hash_leftovers_inserted{$chr}{$pos}{$gene}){
						if(not exists $hash_leftovers{$chr}{$pos}){
			 				$hash_leftovers{$chr}{$pos} = 1;
			 			}
			 			if(exists $hash_leftovers_genes{$chr}{$pos}){
			 				my $array = $hash_leftovers_genes{$chr}{$pos};
			 				push(@$array, $hash_inner);
			 			}else{
			 				my @array = ($hash_inner);
			 				$hash_leftovers_genes{$chr}{$pos} = \@array;
			 			}
			 			$hash_leftovers_inserted{$chr}{$pos}{$gene} = 1;
					}
				}
			}
			
			splice(@$array_genes,0,$quant_remove);
		}
	}
	
	$file_output_leftovers = $file_output.".remaining";
	open(OUT, ">".$file_output_leftovers);
	print OUT sprintf("#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
	"'dbSNP ID'", "Chr", "'Gene Name'","'Gene Start'", "'Gene End'", "'dbSNP pos'", "'QTL db'", "'QTL Category'", "'QTL Start'", "'QTL End'", "'QTL Info'");
	
	my @array_chr = sort {$a<=>$b} keys %hash_leftovers;
	for(@array_chr){ 
		my $chr = $_;
		
		my @array_pos = sort {$a<=>$b} keys %{$hash_leftovers{$chr}};
		
		for(@array_pos){
			my $pos = $_;
			
			my @array_gene_inner = sort{$a->{POS_START}<=>$b->{POS_START}} @{$hash_leftovers_genes{$chr}{$pos}};
			my @closest_markers=();
			my $minimum = 100001;
			for(@array_gene_inner){
				my $hash_gene_inner = $_;
				my $start = $hash_gene_inner->{POS_START};
				my $end = $hash_gene_inner->{POS_END};
				
				my $distance = 0;
				if($pos <= $start ){
					$distance = $start - $pos;
				}else{
					$distance = $pos - $end;
				}
				
				if($distance < $minimum){
					$minimum = $distance;
					@closest_markers=();
					push(@closest_markers, $hash_gene_inner);
				}elsif($distance == $minimum){
					push(@closest_markers, $hash_gene_inner);
				}
			}
			
			for(@closest_markers){
				my $hash_gene_inner = $_;
				my $gene = $hash_gene_inner->{NAME};
				my $start = $hash_gene_inner->{POS_START};
				my $end = $hash_gene_inner->{POS_END};
				print OUT sprintf("%s\t%s\t%s\t%d\t%d\t%d",
									 	 $hash_snp_position_gene{$chr}{$pos}, $chr, $gene,$start, $end, $pos);
				my $hash_qtl_innter = $hash_qtl{$chr}{$hash_snp_position_gene{$chr}{$pos}};
			 	while(my ($key, $value) = each(%$hash_qtl_innter)){
			 		for(@$value){
			 			print OUT sprintf("\t%s\t%s\t%d\t%d\t%s", $_->{QTL}, $key, $_->{POS_START},
			 				$_->{POS_END}, $_->{INFO});
			 		}	
			 	}
				print OUT "\n";
				
				
				if(exists $hash_genes_list{$gene}){
			 		my $category_qtl = "";
			 		my $hash_qtl_innter = $hash_qtl{$chr}{$hash_snp_position_gene{$chr}{$pos}};
				 	while(my ($key, $value) = each(%$hash_qtl_innter)){
				 		if($category_qtl eq ""){
			 				$category_qtl = $key;
			 			}else{
			 				$category_qtl .= ",".$key;
			 			}									 			
				 	}
			 		
					my %hash = (NAME => $gene, CATEGORY_QTL=>$category_qtl, DBSNP_ID=>$hash_snp_position_gene{$chr}{$pos}, FILE=>$file_output_leftovers);
					if(exists $hash_genes_in_files{$gene}){
						my $array = $hash_genes_in_files{$gene};
						push(@$array, \%hash);
					}else{
						my @array=(\%hash);
						$hash_genes_in_files{$gene} = \@array;
					}
				}
				
			}
			
		}
	}
	
	close(OUT);
}

#########################################################################################################
# This function export all genes provided by $file_genes_list if those genes were exported in another file.
#########################################################################################################
sub export_gene_qtl{
	$file_output_gene_qtl = $file_output.".gene_list";
	open(OUT, ">".$file_output_gene_qtl);
	
	print OUT sprintf("#%s\t%s\t%s\t%s\n", "'Gene Name'", "'QTL Category'", "'DBSNP ID'", "'File associated'");
	for(@array_genes_list){
		my $gene = $_;
		if(exists $hash_genes_in_files{$gene}){
			my $array = $hash_genes_in_files{$gene};
			for(@$array){
				print OUT sprintf("%s\t%s\t%s\t%s\n", $_->{NAME}, $_->{CATEGORY_QTL}, $_->{DBSNP_ID}, $_->{FILE});
			}
			
		}
		
	}
	
	close(OUT);
	
	
	
}

#########################################################################################################
# This function reads a dbSNP file and store the information in a hash.
#########################################################################################################
sub read_dbsnp{
	
	open(IN, $file_dbsnp);
	my $count = 0;
	my $wrap_line = 0;
	%hash_count = ();
	while(<IN>){
		chomp;
		unless(/^(\s|\t|#)*$/){
			my @columns = split /\t/;
			my $chr = $columns[0];
			
			if(lc($chr) == 'x'){
				$chr = 30;
			}
			
			my $pos = $columns[1];
			my $gene = $columns[2];
			
			if(exists $hash_qtl{$chr}{$gene}){
				$hash_dbsnp{$chr}{$gene}=$pos;
				$hash_snp_position_gene{$chr}{$pos} = $gene;
			}
		}
		else{
			$count_lines_dbsnp--;
		}
		
		# the next lines are just to calculate and printing the percentage of progress
		my $percent_lines = int(($count++ * 100)/$count_lines_dbsnp);		
		unless(exists $hash_count{$percent_lines}){
			print " | $percent_lines%";
			$wrap_line++;
			$hash_count{$percent_lines} = $percent_lines;
			if( ($percent_lines > 0 and $wrap_line % 10 == 0) or $percent_lines==100 ){
				print "\n";
			}
		}
	}
	unless(exists $hash_count{100}){
		print " | 100%\n";
	}
	
	close(IN);
	
}

#########################################################################################################
# This function validates the parameters provided from user to the script
#########################################################################################################
sub validate_parameters{
	my $allExists = 1;
	my $fileExists = 1;
	
	if ( defined $help ) {
		print usage();
		exit 0;
	}
	
	unless ( defined $file_dbsnp ) {
		$allExists = 0;
	}
	
	unless ( defined $file_gff ) {
		$allExists = 0;
	}
	
	unless ( defined $file_qtl ) {
		$allExists = 0;
	}
	
	unless ( defined $file_output ) {
		$allExists = 0;
	}
	
	if ( defined $file_genes_list ) {
		unless ( -e $file_genes_list ) {
			print STDERR "$file_genes_list doesn't exists.\n";
			$fileExists = 0;
		}
	}
	
	if ($allExists) {
		unless ( -e $file_dbsnp ) {
			print STDERR "$file_dbsnp doesn't exists.\n";
			$fileExists = 0;
		}
		if(-d $file_dbsnp){
			print STDERR "$file_dbsnp is a directory.\n";
			$fileExists = 0;
		}
		
		unless ( -e $file_gff ) {
			print STDERR "$file_gff doesn't exists.\n";
			$fileExists = 0;
		}
		unless ( -e $file_qtl ) {
			print STDERR "$file_qtl doesn't exists.\n";
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
sub usage{
	my $usage = <<FOO;
Usage:
	perl $0 -g gff_file -q qtl_file -d dbsnp_file -o output_file [-h genes_list]
FOO
	return $usage; 
	
}
