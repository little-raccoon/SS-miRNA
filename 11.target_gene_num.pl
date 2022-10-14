#!usr/bin/perl -w
# Guo Zhonglong
# 2021-08-11
# calculate number of target genes

##### INPUT & OUTPUT
my $bg_path="/lustre1/leili_pkuhpc/guozhl/TID/GO/GO_input_for_analysis/";
#my $species_file="/lustre1/leili_pkuhpc/guozhl/TID/species_specific/species_order.lst";
my $species_file="/lustre1/leili_pkuhpc/guozhl/TID/GO/species_for_GO.lst";
my $SS_feature_file="/lustre1/leili_pkuhpc/guozhl/TID/feature/SS_feature.out";
my $NSS_feature_file="/lustre1/leili_pkuhpc/guozhl/TID/feature/NSS_feature.out";
my $target_genes_file="/lustre1/leili_pkuhpc/guozhl/TID/GO/PmiREN2.0_target_genes_add.txt";



my $noTarget_miRNA_out_file="/lustre1/leili_pkuhpc/guozhl/TID/GO/noTarget_miRNA.lst";
##### GLOBAL
my %target_GO_hash;

##### MAIN
### STEP 1: parse target genes file
my %target_hash=parse_targets($target_genes_file);

### STEP 2: extract target genes of SS-miRNA NSS-miRNA
my %SS_miRNA_target_hash=get_targets($SS_feature_file);
my %NSS_miRNA_target_hash=get_targets($NSS_feature_file);

open SPECIES,$species_file or die "$0 fileOpenError: $species_file\n";
while(<SPECIES>){
	chomp;
	my $species=$_;
	if(exists $SS_miRNA_target_hash{$species} && exists $NSS_miRNA_target_hash{$species}){
		my @SS_miRNA_array=keys %{$SS_miRNA_target_hash{$species}};
		my @NSS_miRNA_array=keys %{$NSS_miRNA_target_hash{$species}};
		my $SS_target_total=0;
		my $NSS_target_total=0;
		my $SS_miRNA_total=0;
		my $NSS_miRNA_total=0;
		my $SS_num=@SS_miRNA_array;
		my $NSS_num=@NSS_miRNA_array;
		my $SS_noTarget=0;
		my $NSS_noTarget=0;
		foreach my $miRNA(@SS_miRNA_array){
			if(exists $target_hash{$miRNA}){
				my @targets_array=keys %{$target_hash{$miRNA}};
				my $target_num=@targets_array;
				$SS_target_total+=$target_num;
				$SS_miRNA_total++;
				#print "$species\tSS\t$miRNA\t$target_num\n";
			}
			else{
				#print "SS-miRNA $miRNA no target gene\n";
				$SS_noTarget++;
			}
		}
		foreach my $miRNA(@NSS_miRNA_array){
			if(exists $target_hash{$miRNA}){
				my @targets_array=keys %{$target_hash{$miRNA}};
				my $target_num=@targets_array;
				$NSS_target_total+=$target_num;
				$NSS_miRNA_total++;
				#print "$species\tNSS\t$miRNA\t$target_num\n";
			}
			else{
				#print "NSS-miRNA $miRNA no target gene\n";
				$NSS_noTarget++;
			}
		}
		#print "$species SS $SS_num NSS $NSS_num\n";
		my $SS_ave_targets=$SS_target_total/$SS_miRNA_total;
		my $NSS_ave_targets=$NSS_target_total/$NSS_miRNA_total;
		print "$species\t$SS_ave_targets\t$NSS_ave_targets\n";
		#print "$species\t$SS_num\t$SS_noTarget\t$NSS_num\t$NSS_noTarget\n";
	}
	else{die "$0 species $species not in hash\n";}
	
}
close SPECIES;
#close OUT1;

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
		if($hybrid eq "NA" && $pare eq "NA" && $ps>3){
		
		}
		else{
			my $tmp1=$ps."\t".$hybrid."\t".$pare;
			$hash{$miRNA}{$gene}=$tmp1;
		}
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