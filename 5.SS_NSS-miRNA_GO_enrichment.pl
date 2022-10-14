#!usr/bin/perl -w
# Guo Zhonglong
# 2021-08-05
# extract target genes of SS-miRNA and NSS-miRNA

##### INPUT & OUTPUT
my $bg_path="/lustre1/leili_pkuhpc/guozhl/TID/GO/GO_input_for_analysis/";
my $species_for_GO_file="/lustre1/leili_pkuhpc/guozhl/TID/GO/species_for_GO.lst";
my $SS_feature_file="/lustre1/leili_pkuhpc/guozhl/TID/feature/SS_feature.out";
my $NSS_feature_file="/lustre1/leili_pkuhpc/guozhl/TID/feature/NSS_feature.out";
my $target_genes_file="/lustre1/leili_pkuhpc/guozhl/TID/GO/PmiREN2.0_target_genes_add.txt";



my $noTarget_miRNA_out_file="/lustre1/leili_pkuhpc/guozhl/TID/GO/noTarget_miRNA.lst";
my $SS_targets_out_file="/lustre1/leili_pkuhpc/guozhl/TID/GO/GO_input_for_analysis/SS-miRNA_targets.lst";
my $NSS_targets_out_file="/lustre1/leili_pkuhpc/guozhl/TID/GO/GO_input_for_analysis/NSS-miRNA_targets.lst";
##### GLOBAL
my %target_GO_hash;

##### MAIN
### STEP 1: parse target genes file
my %target_hash=parse_targets($target_genes_file);

### STEP 2: parse GO-target
my @GO_files=glob $bg_path."*out";
foreach my $GO_file(@GO_files){
	open GO_FILE,$GO_file or die "$0 fileOpenError: $GO_file\n";
	while(<GO_FILE>){
		chomp;
		my($gene)=split "\t",$_;
		$target_GO_hash{$gene}=$_;
	}
	close GO_FILE;
}

### STEP 3: extract target genes of SS-miRNA NSS-miRNA
my %SS_miRNA_target_hash=get_targets($SS_feature_file);
my %NSS_miRNA_target_hash=get_targets($NSS_feature_file);

open OUT1,">".$noTarget_miRNA_out_file;
open SS_OUT,">".$SS_targets_out_file;
open NSS_OUT,">".$NSS_targets_out_file;
open SPECIES,$species_for_GO_file or die "$0 fileOpenError: $species_for_GO_file\n";
while(<SPECIES>){
	chomp;
	my $species=$_;
	if(exists $SS_miRNA_target_hash{$species} && exists $NSS_miRNA_target_hash{$species}){
		my @SS_miRNA_array=keys %{$SS_miRNA_target_hash{$species}};
		my @NSS_miRNA_array=keys %{$NSS_miRNA_target_hash{$species}};
		foreach my $miRNA(@SS_miRNA_array){
			if(exists $target_hash{$miRNA}){
				my @targets_array=keys %{$target_hash{$miRNA}};
				#$SS_num++;
				foreach my $gene(@targets_array){
					if(exists $target_GO_hash{$gene}){
						print SS_OUT "$gene\n";
					}
					else{
						#print "SS-miRNA $species $miRNA $gene no GO term\n";
					}
				}
			}
			else{
				print OUT1 "SS-miRNA $miRNA no target gene\n";
			}
		}
		foreach my $miRNA(@NSS_miRNA_array){
			if(exists $target_hash{$miRNA}){
				my @targets_array=keys %{$target_hash{$miRNA}};
				#$NSS_num++;
				foreach my $gene(@targets_array){
					if(exists $target_GO_hash{$gene}){
						print NSS_OUT "$gene\n";
					}
					else{
						#print "NSS-miRNA $species $miRNA $gene no GO term\n";
					}
				}
			}
			else{
				print OUT1 "NSS-miRNA $miRNA no target gene\n";
			}
		}
		#print "$species SS $SS_num NSS $NSS_num\n";
	}
	else{die "$0 species $species not in hash\n";}
	
}
close SPECIES;
close OUT1;
close SS_OUT;
close NSS_OUT;

##### FUNCTION
sub parse_targets{
	my %hash;
	my($infile)=@_;
	open INFILE,$infile or die "$0 fileOpenError: $infile\n";
	while(<INFILE>){
		chomp;
		my @line=split "\t",$_;
		my($miRNA,$gene,$ps,$hybrid,$pare)=@line[1..5];
		$miRNA=~s/miR/MIR/g;
		my $tmp1=$ps."\t".$hybrid."\t".$pare;
		$hash{$miRNA}{$gene}=$tmp1;
	}
	close INFILE;
	return %hash;
}

sub get_targets{
	my %hash;
	my($infile)=@_;
	open INFILE,$infile or die "$0 fileOpenError: $infile\n";
	while(<INFILE>){
		chomp;
		my($species,$miRNA)=split "\t",$_;
		$hash{$species}{$miRNA}=1;
	}
	close INFILE;
	return %hash;
}