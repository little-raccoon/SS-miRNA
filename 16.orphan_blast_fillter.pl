#!/usr/bin/perl

my $Osa_pre_file="/lustre1/leili_pkuhpc/guozhl/TID/orphan/Osa_miRNA_pre.fa";
my $blastout_path="/lustre1/leili_pkuhpc/guozhl/TID/orphan/blastout/";

open OSA,$Osa_pre_file or die "$0 fileOpenError: $Osa_pre_file\n";
while(<OSA>){
	chomp;
	my $line=$_;
	if($line=~/^>/){
		my $miRNA=$line;
		$miRNA=~s/>//g;
		my $seq=<OSA>;
		chomp $seq;
		my $len=length($seq);
		$miRNA_len_hash{$miRNA}=$len;
		#print "$miRNA\t$len\n";
	}
}
close OSA;

my @blastout_files=glob $blastout_path."*blastout";
foreach my $file(@blastout_files){
	#print "$file\n";
	my $threshold=90; #匹配的阈值
	my $threshold_ratio=$threshold/100;
	my $outfile=$file.".$threshold.filterout";
	open INFILE,$file or die "fileOpenError: $file\n";
	open OUT,">".$outfile;
	while(<INFILE>){
		chomp;
		my $line=$_;
		my @line=split "\t",$line;
		my $miRNA=$line[0];
		my $match_len=$line[3];
		#print "$miRNA\t$match_len\n";
		if(exists $miRNA_len_hash{$miRNA}){
			my $miRNA_len=$miRNA_len_hash{$miRNA};
			my $ratio=$match_len/$miRNA_len;
			if($ratio>=$threshold_ratio){
				#print "$ratio\t$threshold_ratio\n";
				print OUT "$ratio\t$line\n";
			}
		}
		else{
			die "$0: $miRNA not in miRNA_len_hash\n";
		}
	}
	close INFILE;
	close OUT;
}
