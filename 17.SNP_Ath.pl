#!/usr/bin/perl

my $VCF_file="/gpfs1/leili_pkuhpc/guozhl/SNP/Ath1001/1001genomes_snp-short-indel_with_tair10_only_ACGTN.vcf";
my $basic_file="/lustre1/leili_pkuhpc/guozhl/TID/PmiREN2.0_data/PmiREN2.0_basic_info.txt";
my $out_path="/lustre1/leili_pkuhpc/guozhl/TID/SNP/Ath/";

open BASIC,$basic_file or die "$0 fileOpenError: $basic_file\n";
while(<BASIC>){
	chomp;
	my $line=$_;
	my @line=split "\t",$line;
	my($species,$miRNA_fam,$chr,$start,$end)=@line[1..5];
	my $miRNA=$line[9];
	#print "$species\t$miRNA_fam\t$chr\t$miRNA\n";
	if($species eq "Arabidopsis thaliana"){
		#print "$species\n";
		$chr=~s/Ath-Chr//g;
		#print "$chr\n";
		my $out=$out_path.$miRNA.".vcfout";
		system "vcftools --vcf $VCF_file --chr $chr --from-bp $start --to-bp $end --out $out --counts";
	}
}
close BASIC;