#!usr/bin/perl -w
# Guo Zhonglong
# 2021-08-10
# extract expression value of SS-miRNAs and NSS-miRNAs

##### INPUT & OUTPUT
my $expression_file="/lustre1/leili_pkuhpc/guozhl/TID/expression/wuzhongRNAMatureMap.txt";
my $SS_file="/lustre1/leili_pkuhpc/guozhl/TID/feature/SS_feature.out";
my $NSS_file="/lustre1/leili_pkuhpc/guozhl/TID/feature/NSS_feature.out";

##### MAIN
### STEP 1：parse expression file
my %exp_hash=parse_expression($expression_file);

### STEP 2: extract specific species
foreach my $miRNA(sort keys %exp_hash){
	my @tissues=sort keys %{$exp_hash{$miRNA}};
	my $tissue_num=@tissues;
	print "$miRNA\t$tissue_num\t@tissues\n";
}

###### FUNCTION
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
		$exp=~s/},/}/g;
		#$exp=~s/://g;
		#$exp=~s/,//g;
		#$exp=~s/[a-zA-Z]*//g;
		$exp=~s/\{//g;
		my @exp=split "}",$exp;
		#print "$miRNA";
		foreach my $term(@exp){
			if($term=~/name:(.*),value:(.*)/){
				my $tissue=$1;
				my $value=$2;
				$hash{$miRNA}{$tissue}=$value;
				#print "$tissue-$value";
			}
			else{die "$0: reg express error in $miRNA $term\n";}
		}
		#print "$miRNA\t$ave_rpm\n";
	}
	close INFILE;
	return %hash;
}
