#!/usr/bin/perl
# Guo Zhonglong
# 2021-09-02
# find different types of nt pairs

##### INPUT & OUTPUT
my $fold_file="/lustre1/leili_pkuhpc/guozhl/TID/nt_pair/PmiREN2_pre.rnafold";

##### MAIN

### STEP 1: parse fold file
my %fold_hash=parse_fold($fold_file);

foreach my $miRNA(sort keys %fold_hash){
	my $tmp1=$fold_hash{$miRNA};
	my($seq,$struc)=split ";",$tmp1;
}

##### FUNCTION
sub parse_fold{
	my($infile)=@_;
	open INFILE,$infile or die "$0 fileOpenError: $infile\n";
	my %hash;
	while(<INFILE>){
		chomp;
		my $line=$_;
		if($line=~/^>/){
			my $miRNA=$line;
			$miRNA=~s/>//g;
			my $seq=uc <INFILE>;
			my $struc=<INFILE>;
			chomp $seq;
			chomp $struc;
			$seq=~s/T/U/g;
			$struc=~s/ .*//g;
			#print "$struc\n";
			$hash{$miRNA}=$seq.";".$struc;
		}
	}
	close INFILE;
	return %hash;
}
