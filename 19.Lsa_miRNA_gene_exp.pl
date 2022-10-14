#!/usr/bin/perl

my $miRNA_exp_file=shift;
my $gene_exp_file=shift;
my $miRNA_lst_file=shift;

open MIRNA,$miRNA_exp_file;
while(<MIRNA>){
	chomp;
	my $line=$_;
	my @line=split "\t",$line;
	$miRNA_hash{$line[0]}=$line;
}
close MIRNA;

open GENE,$gene_exp_file;
while(<GENE>){
	chomp;
	my $line=$_;
	my @line=split "\t",$line;
	my $gene=$line[0];
	$gene=~s/\..*//g;
	$gene_hash{$gene}{$line[0]}=$line;
}
close GENE;

open LST,$miRNA_lst_file;
while(<LST>){
	chomp;
	my $line=$_;
	my ($miRNA,$gene)=split "\t",$line;
	if(exists $miRNA_hash{$miRNA} && exists $gene_hash{$gene}){
		my @tmp=sort keys %{$gene_hash{$gene}};
		#print "@tmp\n";
		#print "$out\n";
		foreach my $transcript(@tmp){
			my $out=$miRNA_hash{$miRNA}."\t".$gene."\t".$gene_hash{$gene}{$transcript};
			print "$out\n";
		}
	}
	else{
		#print "$line\n";
	}
}
close LST;