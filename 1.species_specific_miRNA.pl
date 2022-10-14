#!usr/bin/perl -w
# Guo Zhonglong
# 2021-07-29
# identify species-specific miRNAs using data in PmiREN2.0

##### INPUT & OUTPUT
my $PmiREN2_species_file="/lustre1/leili_pkuhpc/guozhl/TID/PmiREN2.0_data/wuzhong.txt";
my $PmiREN2_basic_file="/lustre1/leili_pkuhpc/guozhl/TID/PmiREN2.0_data/wuzhongRNA.txt";
my $sam_file="/lustre1/leili_pkuhpc/guozhl/TID/species_specific/PmiREN2_mature.sam";

my $SS_miRNA_out_file="/lustre1/leili_pkuhpc/guozhl/TID/species_specific/SS-miRNA_basic_info.out";
my $NSS_miRNA_out_file="/lustre1/leili_pkuhpc/guozhl/TID/species_specific/NSS-miRNA_basic_info.out";
my $SS_percent_out_file="/lustre1/leili_pkuhpc/guozhl/TID/species_specific/SS-miRNA_percent.out";

##### GLOBAL
my %miRNA_spe_hash;
my %miRNA_basic_hash;
my %spe_miRNA_num_hash;
my %NSS_miRNA_hash; # NSS-miRNA means non-species-specific miRNA
my %SS_miRNA_hash; # SS-miRNA means species-specific miRNA




##### MAIN
### STEP 1: parse data in PmiREN2.0
my($miRNA_spe_hash,$miRNA_basic_hash,$spe_miRNA_num_hash)=parse_PmiREN($PmiREN2_species_file,$PmiREN2_basic_file);
%miRNA_spe_hash=%{$miRNA_spe_hash};
%miRNA_basic_hash=%{$miRNA_basic_hash};
%spe_miRNA_num_hash=%{$spe_miRNA_num_hash};

### STEP 2: parse sam file and find NSS-miRNA and SS-miRNA
%NSS_miRNA_hash=parse_sam($sam_file,\%miRNA_spe_hash);
foreach my $miRNA(keys %miRNA_spe_hash){
	if(! exists $NSS_miRNA_hash{$miRNA}){
		$SS_miRNA_hash{$miRNA}=$miRNA_spe_hash{$miRNA};
	}
}

### STEP 3: out basic info of SS-miRNA and NSS-miRNA
open SS_MIRNA_OUT,">".$SS_miRNA_out_file;
foreach my $SS_miRNA(sort keys %SS_miRNA_hash){
	if(exists $miRNA_basic_hash{$SS_miRNA}){
		#print "$miRNA_basic_hash{$SS_miRNA}\n";
		print SS_MIRNA_OUT "$miRNA_basic_hash{$SS_miRNA}\n";
	}
	else{die "$0 SS_miRNA: $SS_miRNA not in hash: miRNA_basic_hash\n";}
}
close SS_MIRNA_OUT;

open NSS_MIRNA_OUT,">".$NSS_miRNA_out_file;
foreach my $NSS_miRNA(sort keys %NSS_miRNA_hash){
	if(exists $miRNA_basic_hash{$NSS_miRNA}){
		print NSS_MIRNA_OUT "$miRNA_basic_hash{$NSS_miRNA}\n";
	}
	else{die "$0 NSS_miRNA: $NSS_miRNA not in hash: miRNA_basic_hash\n";}
}
close NSS_MIRNA_OUT;


### STEP 4: calculate percent
foreach my $SS_miRNA(keys %SS_miRNA_hash){
	my $species=$SS_miRNA_hash{$SS_miRNA};
	$spe_SSmiRNA_num_hash{$species}++;
}

open SS_PERCENT_OUT,">".$SS_percent_out_file;
foreach my $species(sort keys %spe_miRNA_num_hash){
	my $SS_percent;
	my $total_miRNA_num=$spe_miRNA_num_hash{$species};
	#print "$species\t$total_miRNA_num\n";
	if(exists $spe_SSmiRNA_num_hash{$species}){
		$SS_percent=$spe_SSmiRNA_num_hash{$species}/$total_miRNA_num*100;
	}
	else{
		$SS_percent=0;
	}
	print SS_PERCENT_OUT "$species\t$total_miRNA_num\t$SS_percent\n";
}
close SS_PERCENT_OUT;





##### FUNCTION
sub parse_PmiREN{
	my($PmiREN2_species_file,$PmiREN2_basic_file)=@_;
	my %spe_hash;
	my %miRNA_spe_hash;
	open SPE,$PmiREN2_species_file or die "$0 fileOpenError: $PmiREN2_species_file\n";
	while(<SPE>){
		chomp;
		my @line=split "\t",$_;
		my($id,$species)=($line[0],$line[2]);
		$spe_hash{$id}=$species;
	}
	close SPE;
	open BASIC,$PmiREN2_basic_file or die "$0 fileOpenError: $PmiREN2_basic_file\n";
	while(<BASIC>){
		chomp;
		my $species;
		my $line=$_;
		my @line=split "\t",$line;
		my($miRNA,$family,$mature,$species_id)=($line[13],$line[2],$line[18],$line[1]);
		if(exists $spe_hash{$species_id}){
			$species=$spe_hash{$species_id};
		}
		else{die "$0 $species_id not in hash: spe_hash\n"}
		$miRNA_spe_hash{$miRNA}=$species;
		$miRNA_basic_hash{$miRNA}=$line;
		#print "$line\n";
		$spe_miRNA_num_hash{$species}++;
		#print "$species\n";
	}
	close BASIC;
	return \%miRNA_spe_hash,\%miRNA_basic_hash,\%spe_miRNA_num_hash;
}

sub parse_sam{
	my %NSS_miRNA_hash;
	my($sam_file,$hash)=@_;
	open SAM,$sam_file or die "$0 fileOpenError: $sam_file\n";
	while(<SAM>){
		chomp;
		if($_!~/^@/){
			my @line=split "\t",$_;
			my($query,$strand,$db,$site)=@line[0..3];
			if(exists $$hash{$query} && exists $$hash{$db}){
				my $query_spe=$$hash{$query};
				my $db_spe=$$hash{$db};
				if($query_spe ne $db_spe){
					#print "$query_spe\t$query\t$db_spe\t$db\t$strand\t$site\n";
					$NSS_miRNA_hash{$query}=$query_spe;
					$NSS_miRNA_hash{$db}=$db_spe;
				 }
			}
			else{die "$0 $query or $db not in hash: spe_hash\n"}
		}
	}
	close SAM;
	return %NSS_miRNA_hash;
}