#!usr/bin/perl -w
# Guo Zhonglong
# 2021-08-06
# prepare the input file for machine learning

##### INPUT & OUTPUT

my $SS_feature_file="/lustre1/leili_pkuhpc/guozhl/TID/ML/SS_feature.out";
my $NSS_feature_file="/lustre1/leili_pkuhpc/guozhl/TID/ML/NSS_feature.out";

my $SS_NSS_feature_for_ML_out_file="/lustre1/leili_pkuhpc/guozhl/TID/ML/SS_NSS_feature_for_ML.out";
##### GLOBAL

##### MAIN
### STEP 1: parse origin file
my %SS_feature_hash=parse_orgin($SS_feature_file,"1"); # type 1=>SS-miRNA; type 0=>NSS-miRNA
my %NSS_feature_hash=parse_orgin($NSS_feature_file,"0");
### STEP 2: merge SS-miRNA and NSS-miRNA features
open OUT,">".$SS_NSS_feature_for_ML_out_file;
print OUT "species\tmiRNA\ttype\tmature_len\tpre_len\tfirstBase\tm_A\tm_T\tm_C\tm_G\tp_A\tp_T\tp_C\tp_G\tnmef\n";
foreach my $k(sort keys %SS_feature_hash){
	print OUT "$SS_feature_hash{$k}\n";
}
foreach my $k(sort keys %NSS_feature_hash){
	print OUT "$NSS_feature_hash{$k}\n";
}
close OUT;
##### FUNCTION
sub parse_orgin{
	my %hash;
	my($infile,$type)=@_;
	open INFILE,$infile or die "$0 fileOpenError: $infile\n";
	while(<INFILE>){
		chomp;
		my $mark;
		my $line=$_;
		my @line=split "\t",$line;
		my $firstBase=$line[4];
		if($firstBase eq "A" || $firstBase eq "a"){$mark=1;}
		elsif($firstBase eq "T" || $firstBase eq "t" || $firstBase eq "U" || $firstBase eq "t"){$mark=2;}
		elsif($firstBase eq "C" || $firstBase eq "c"){$mark=3;}
		elsif($firstBase eq "G" || $firstBase eq "g"){$mark=4;}
		else{die "$line first base error\n";}
		$line[4]=$mark;
		$line[5]=$line[5]/$line[2];
		$line[6]=$line[6]/$line[2];
		$line[7]=$line[7]/$line[2];
		$line[8]=$line[8]/$line[2];
		$line[9]=$line[9]/$line[3];
		$line[10]=$line[10]/$line[3];
		$line[11]=$line[11]/$line[3];
		$line[12]=$line[12]/$line[3];
		my $tmp1=join "\t",@line[0..1];
		my $tmp2=join "\t",@line[2..$#line];
		$hash{$line[1]}=$tmp1."\t".$type."\t".$tmp2;
		#print "$tmp1\t$type\t$tmp2\n";
	}
	close INFILE;
	return %hash;
}