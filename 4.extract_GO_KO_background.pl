#!usr/bin/perl -w
# Guo Zhonglong
# 2021-08-05
# extract GO background

##### INPUT & OUTPUT
my $phytozome_files_path="/lustre1/leili_pkuhpc/guozhl/TID/GO/phytozome_files_path.txt";

my $GO_origin_out_path="/lustre1/leili_pkuhpc/guozhl/TID/GO/GO_origin_out/";
my $GO_input_for_analysis="/lustre1/leili_pkuhpc/guozhl/TID/GO/GO_input_for_analysis/";

##### GLOBAL


##### MAIN
### STEP 1: parse GO annotation files
open FILES,$phytozome_files_path or die "$0 fileOpenError: $phytozome_files_path\n";
while(<FILES>){
	chomp;
	if($.>1){
		my($species,$genome_file,$mRNA_file,$cds_file,$protein_file,$GO_file,$gff3_file)=split "\t",$_;
		#print "$species\t$GO_file\n";
		$species=~s/ /_/g;
		my $GO_origin_out_file=$GO_origin_out_path.$species.".origin.out";
		my $GO_background_out_file=$GO_input_for_analysis.$species.".GO.background.out";
		my $KO_background_out_file=$GO_input_for_analysis.$species.".KO.background.out";
		open GO_OUT1,">".$GO_origin_out_file;
		open GO_OUT2,">".$GO_background_out_file;
		open KO_OUT1,">".$KO_background_out_file;
		open GO,$GO_file or die "$0 fileOpenError: $GO_file\n";
		while(<GO>){
			chomp;
			if($.>1){
				my @line=split "\t",$_;
				my($gene,$ko,$go,$symbol,$define)=("na","na","na","na","na");
				$gene=$line[2];
				if($line[8]){$ko=$line[8];}
				if($line[9]){$go=$line[9];}
				if($line[11]){$symbol=$line[11];}
				if($line[12]){$define=$line[12];}
				#($gene,$ko,$go,$symbol,$define)=($line[2],$line[8],$line[9],$line[11],$line[12]);
				print GO_OUT1 "$gene\t$ko\t$go\t$symbol\t$define\n";
				my @ko=split ",",$ko;
				my @go=split ",",$go;
				foreach my $k(@ko){print KO_OUT1 "$gene\t$k\n";}
				foreach my $k(@go){print GO_OUT2 "$gene\t$k\n";}
			}
		}
		close GO;
		close GO_OUT1;
		close GO_OUT2;
		close KO_OUT1;
	}
}
close FILES;

##### FUNCTION
