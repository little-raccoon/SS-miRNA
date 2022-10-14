#!usr/bin/perl -w
# Guo Zhonglong
# 2021-07-30
# group SS-miRNA into family

##### INPUT & OUTPUT
my $PmiREN2_species_file="/lustre1/leili_pkuhpc/guozhl/TID/PmiREN2.0_data/wuzhong.txt";
my $PmiREN2_basic_file="/lustre1/leili_pkuhpc/guozhl/TID/PmiREN2.0_data/wuzhongRNA.txt";
my $SS_basic_file="/lustre1/leili_pkuhpc/guozhl/TID/species_specific/SS-miRNA_basic_info.out";


##### GLOBAL



##### MAIN
### STEP 1: parse SS_basic_file
my %All_fam_hash=parse_basic($PmiREN2_species_file,$PmiREN2_basic_file);
my %SS_fam_hash=parse_basic($PmiREN2_species_file,$SS_basic_file);

### STEP 2: calculate number of miRNA families in individual species
foreach my $species(sort keys %All_fam_hash){
	my @all_fam_array=keys %{$All_fam_hash{$species}};
	my $all_fam_num=@all_fam_array;
	my @SS_fam_array;
	my $SS_fam_num=0;
	if(exists $SS_fam_hash{$species}){
		@SS_fam_array=keys %{$SS_fam_hash{$species}};
		$SS_fam_num=@SS_fam_array;
	}
	my $SS_fam_percent=$SS_fam_num/$all_fam_num*100;
	print "$species\t$SS_fam_num\t$all_fam_num\t$SS_fam_percent\t@SS_fam_array\n";
}

### STEP 3:



##### FUNCTION
sub parse_basic{
	my($PmiREN2_species_file,$basic_file)=@_;
	my %spe_hash;
	my %spe_fam_hash;
	open SPE,$PmiREN2_species_file or die "$0 fileOpenError: $PmiREN2_species_file\n";
	while(<SPE>){
		chomp;
		my @line=split "\t",$_;
		my($id,$species)=($line[0],$line[2]);
		$spe_hash{$id}=$species;
	}
	close SPE;
	open BASIC,$basic_file or die "$0 fileOpenError: $basic_file\n";
	while(<BASIC>){
		chomp;
		my $species;
		my $line=$_;
		my @line=split "\t",$line;
		my($miRNA,$family,$mature,$species_id)=($line[13],$line[2],$line[18],$line[1]);
		if(exists $spe_hash{$species_id}){
			$species=$spe_hash{$species_id};
		}
		else{die "$0 $species_id not in hash: spe_hash\n";}
		$spe_fam_hash{$species}{$family}++; # members count in a family in a species
	}
	close BASIC;
	return %spe_fam_hash;
}