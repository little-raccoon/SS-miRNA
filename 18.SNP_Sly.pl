#!/usr/bin/perl

my $VCF_file="/gpfs1/leili_pkuhpc/guozhl/SNP/Sly360/Sly_miRNA_by_li/360_merged_2.50.vcf";
my $basic_file="/lustre1/leili_pkuhpc/guozhl/TID/SNP/Sly/Sly_2.5.bed";
my $out_path="/lustre1/leili_pkuhpc/guozhl/TID/SNP/Sly/";

open BASIC,$basic_file or die "$0 fileOpenError: $basic_file\n";
while(<BASIC>){
	chomp;
	my $line=$_;
	my @line=split "\t",$line;
	my($chr,$start,$end,$miRNA,$strand)=@line;
	#print "$species\t$miRNA_fam\t$chr\t$miRNA\n";
	#print "$species\n";
	#print "$chr\n";
	my $out=$out_path.$miRNA.".vcfout";
	system "vcftools --vcf $VCF_file --chr $chr --from-bp $start --to-bp $end --out $out --counts";
}
close BASIC;