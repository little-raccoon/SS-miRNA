#!usr/bin/perl -w
# Guo Zhonglong
# 2021-08-08
# group SS-miRNA into family

##### INPUT & OUTPUT
my $PmiREN2_species_file="/lustre1/leili_pkuhpc/guozhl/TID/PmiREN2.0_data/wuzhong.txt";
my $SS_basic_file="/lustre1/leili_pkuhpc/guozhl/TID/species_specific/SS-miRNA_basic_info.out";
my $NSS_basic_file="/lustre1/leili_pkuhpc/guozhl/TID/species_specific/NSS-miRNA_basic_info.out";

my $outpath="/lustre1/leili_pkuhpc/guozhl/TID/location/miRNA_bed/";
##### GLOBAL



##### MAIN
### STEP 1: convert basic file to bed file
my($SS_mature_bed_hash,$SS_pre_bed_hash)=parse_basic($PmiREN2_species_file,$SS_basic_file);
my($NSS_mature_bed_hash,$NSS_pre_bed_hash)=parse_basic($PmiREN2_species_file,$NSS_basic_file);
my %SS_mature_bed_hash=%{$SS_mature_bed_hash};
my %SS_pre_bed_hash=%{$SS_pre_bed_hash};
my %NSS_mature_bed_hash=%{$NSS_mature_bed_hash};
my %NSS_pre_bed_hash=%{$NSS_pre_bed_hash};

### STEP 2: output
output(\%SS_mature_bed_hash,"SS_mature",$outpath);
output(\%SS_pre_bed_hash,"SS_pre",$outpath);
output(\%NSS_mature_bed_hash,"NSS_mature",$outpath);
output(\%NSS_pre_bed_hash,"NSS_pre",$outpath);

##### FUNCTION
sub parse_basic{
	my($PmiREN2_species_file,$basic_file)=@_;
	my %spe_hash;
	my %mature_bed_hash;
	my %pre_bed_hash;
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
		my %spe_position_hash;
		my $line=$_;
		my @line=split "\t",$line;
		my($species_id,$chr,$miRNA,$mature_S,$mature_E,$pre_S,$pre_E,$strand)=($line[1],$line[3],$line[13],$line[16],$line[17],$line[11],$line[12],$line[6]);
		if($chr ne "n\.a\."){
			if(exists $spe_hash{$species_id}){
				$species=$spe_hash{$species_id};
			}
			else{die "$0 $species_id not in hash: spe_hash\n";}
			$chr=~s/^...-//g;
			my $mature_bed=$chr."\t".$mature_S."\t".$mature_E."\t".$miRNA."\t.\t".$strand;
			my $pre_bed=$chr."\t".$pre_S."\t".$pre_E."\t".$miRNA."\t.\t".$strand;
			$mature_bed_hash{$species}{$miRNA}=$mature_bed;
			$pre_bed_hash{$species}{$miRNA}=$pre_bed;
		}
	}
	close BASIC;
	return \%mature_bed_hash,\%pre_bed_hash;
}

sub output{
	my($hash,$type,$path)=@_;
	foreach my $species(sort keys %{$hash}){
		my $outfile=$path.$species."_".$type.".bed";
		$outfile=~s/ /_/g;
		open OUT,">".$outfile;
		my @miRNAs=sort keys %{$$hash{$species}};
		foreach my $miRNA(@miRNAs){
			print OUT "${$$hash{$species}}{$miRNA}\n";
		}
		close OUT;
	}
}