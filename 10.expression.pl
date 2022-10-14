#!usr/bin/perl -w
# Guo Zhonglong
# 2021-08-10
# extract expression value of SS-miRNAs and NSS-miRNAs

##### INPUT & OUTPUT
my $expression_file="/lustre1/leili_pkuhpc/guozhl/TID/expression/wuzhongRNAMatureMap.txt";
my $SS_file="/lustre1/leili_pkuhpc/guozhl/TID/feature/SS_feature.out";
my $NSS_file="/lustre1/leili_pkuhpc/guozhl/TID/feature/NSS_feature.out";

##### MAIN

### STEP 1: parse expression file
my %expression_hash=parse_expression($expression_file);

### STEP 2: extract expression values of SS-miRNAs and NSS-miRNAs
my %SS_exp_hash=extract_exp($SS_file,\%expression_hash,"SS");
my %NSS_exp_hash=extract_exp($NSS_file,\%expression_hash,"NSS");

### STEP 3: output
foreach my $species(sort keys %SS_exp_hash){
	print "$species\t$SS_exp_hash{$species}\t$NSS_exp_hash{$species}\n";
}
##### FUNCTION
sub parse_expression{
	my($infile)=@_;
	my %hash;
	open INFILE,$infile or die "$0 fileOpenError: $infile\n";
	while(<INFILE>){
		chomp;
		my @line=split "\t",$_;
		my $miRNA=$line[1];
		my $exp=$line[2];
		$exp=~s/"//g;
		$exp=~s/\[//g;
		$exp=~s/\]//g;
		$exp=~s/://g;
		$exp=~s/,//g;
		$exp=~s/[a-zA-Z]*//g;
		$exp=~s/\{//g;
		my @exp=split "}",$exp;
		my $total_value=0;
		foreach my $value(@exp){
			$total_value+=$value;
		}
		my $count=@exp;
		my $ave_rpm=$total_value/$count;
		$hash{$miRNA}=$ave_rpm;
		#print "$miRNA\t$ave_rpm\n";
	}
	close INFILE;
	return %hash;
}

sub extract_exp{
	my($infile,$hash,$type)=@_;
	my %out_hash;
	my %exp_hash;
	my %count_hash;
	open INFILE,$infile or die "$0 fileOpenError: $infile\n";
	while(<INFILE>){
		chomp;
		my @line=split "\t",$_;
		my($species,$miRNA)=@line;
		#print "$species\t$miRNA\n";
		if(exists ${$hash}{$miRNA}){
			my $ave_rmp=${$hash}{$miRNA};
			#print "$miRNA\t$ave_rmp\t$type\t$species\n";
			$exp_hash{$species}+=$ave_rmp;
			$count_hash{$species}++;
		}
		else{
			#die "$0 species: $species miRNA: $miRNA not in expression hash\n";
		}
	}
	close INFILE;
	foreach my $species(sort keys %exp_hash){
		my $exp=$exp_hash{$species};
		my $count=$count_hash{$species};
		my $spe_ave_exp=$exp/$count;
		$out_hash{$species}=$spe_ave_exp;
	}
	return %out_hash;
}