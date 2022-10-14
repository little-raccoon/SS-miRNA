#!usr/bin/perl -w
# Guo Zhonglong
# 2021-08-10
# find the closest gene of miRNA

##### INPUT & OUTPUT
my $name_file="/lustre1/leili_pkuhpc/guozhl/TID/location/phytozome_files_path.txt";

##### MAIN

### STEP 1: match miRNA bed and gene bed
open FILES,$name_file or die "$0 fileOpenError: $name_file\n";
while(<FILES>){
	chomp;
	my @line=split "\t",$_;
	if($.>1){
		my $species=$line[0];
		my $gff_file=$line[-1];
		$species=~s/ /_/g;
		my $gene_gff_file;
		if($gff_file=~/([^\/]*\.)gff3/){
			$gene_gff_file="/lustre1/leili_pkuhpc/guozhl/TID/location/gene_bed/".$1."bed";
			#print "$gene_gff_file\n";
		}
		else{
			die "$0 no match reg express\n";
		}
		my $SS_gff_file="/lustre1/leili_pkuhpc/guozhl/TID/location/miRNA_bed/".$species."_SS_mature.bed";
		my $NSS_gff_file="/lustre1/leili_pkuhpc/guozhl/TID/location/miRNA_bed/".$species."_NSS_mature.bed";
		if(-e $gene_gff_file && -e $SS_gff_file && -e $NSS_gff_file){
			my $SS_out_file="/lustre1/leili_pkuhpc/guozhl/TID/location/closest_out/".$species.".SS.closest.out";
			my $NSS_out_file="/lustre1/leili_pkuhpc/guozhl/TID/location/closest_out/".$species.".NSS.closest.out";
			`bedtools closest -D b -a $SS_gff_file -b $gene_gff_file > $SS_out_file`;
			`bedtools closest -D b -a $NSS_gff_file -b $gene_gff_file > $NSS_out_file`;
		}
		else{die "$0 species: $species a file not exists\n";}
	}
	
}
close FILES;

##### FUNCTION