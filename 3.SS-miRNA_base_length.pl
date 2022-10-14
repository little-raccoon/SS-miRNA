#!usr/bin/perl -w
# Guo Zhonglong
# 2021-08-01
# calculate features of SS-miRNA and NSS-miRNA, such as base composition, length etc.

##### INPUT & OUTPUT
my $PmiREN2_species_file="/lustre1/leili_pkuhpc/guozhl/TID/PmiREN2.0_data/wuzhong.txt";
my $PmiREN2_basic_file="/lustre1/leili_pkuhpc/guozhl/TID/PmiREN2.0_data/wuzhongRNA.txt";
my $SS_basic_file="/lustre1/leili_pkuhpc/guozhl/TID/species_specific/SS-miRNA_basic_info.out";
my $species_order_file="/lustre1/leili_pkuhpc/guozhl/TID/species_specific/species_order.lst";

my $SS_feature_out_file="/lustre1/leili_pkuhpc/guozhl/TID/feature/SS_feature.out";
my $NSS_feature_out_file="/lustre1/leili_pkuhpc/guozhl/TID/feature/NSS_feature.out";
my $feature_out_file="/lustre1/leili_pkuhpc/guozhl/TID/feature/SS_NSS_features_average.out";

##### GLOBAL
my %species_hash;
my %feature_hash;
my %SS_feature_hash;
my %NSS_feature_hash;


##### MAIN
### STEP 1: parse SS_basic_file
my %All_miRNA_hash=parse_basic($PmiREN2_species_file,$PmiREN2_basic_file);
my %SS_miRNA_hash=parse_basic($PmiREN2_species_file,$SS_basic_file);

### STEP 2: calculate features and store in hash
foreach my $species(sort keys %All_miRNA_hash){
	foreach my $miRNA(sort keys %{$All_miRNA_hash{$species}}){
		my($mature,$pre)=split "\t",$All_miRNA_hash{$species}{$miRNA};
		#print "$species\t$miRNA\t$mature\t$pre\n";
		my $mature_len=length($mature);
		my $pre_len=length($pre);
		my $firstBase=substr($mature,0,1);
		my($A,$T,$C,$G)=base_composite($mature);
		my($pre_A,$pre_T,$pre_C,$pre_G)=base_composite($pre);
		#print "$miRNA\t$A\t$T\t$C\t$G\n";
		my $tmp1=$mature_len."\t".$pre_len."\t".$firstBase;
		my $tmp2=$A."\t".$T."\t".$C."\t".$G;
		my $tmp3=$pre_A."\t".$pre_T."\t".$pre_C."\t".$pre_G;
		$feature_hash{$species}{$miRNA}=$tmp1."\t".$tmp2."\t".$tmp3;
	}
}

### STEP 3: extract features of SS-miRNA and NSS-miRNA
open SS_OUT,">".$SS_feature_out_file;
open NSS_OUT,">".$NSS_feature_out_file;
open SPECIES,$species_order_file or die "$0 fileOpenError: $species_order_file\n";
while(<SPECIES>){
	chomp;
	my $species=$_;
	if(exists $feature_hash{$species}){
		foreach my $miRNA(sort keys %{$feature_hash{$species}}){
			my $tmp1=$feature_hash{$species}{$miRNA};
			if(exists $SS_miRNA_hash{$species}{$miRNA}){ # whether or not a SS-miRNA
				print SS_OUT "$species\t$miRNA\t$tmp1\n";
			}
			else{
				print NSS_OUT "$species\t$miRNA\t$tmp1\n";
			}
		}
	}
	else{die "$0 species $species not in hash: feature_hash\n";}
}
close SPECIES;
close SS_OUT;
close NSS_OUT;

### STEP 4: calculate features in a species
%SS_feature_hash=cal_feature_inSpecies($SS_feature_out_file);
%NSS_feature_hash=cal_feature_inSpecies($NSS_feature_out_file);

open FEATURE_OUT,">".$feature_out_file;
print FEATURE_OUT "species\tSS_ave_mature_len\tSS_ave_pre_len\tSS_ave_A\tSS_ave_T\tSS_ave_C\tSS_ave_G\tSS_ave_pre_A\tSS_ave_pre_T\tSS_ave_pre_C\tSS_ave_pre_G\tSS_ave_firstBase_A\tSS_ave_firstBase_T\tSS_ave_firstBase_C\tSS_ave_firstBase_G\tNSS_ave_mature_len\tNSS_ave_pre_len\tNSS_ave_A\tNSS_ave_T\tNSS_ave_C\tNSS_ave_G\tNSS_ave_pre_A\tNSS_ave_pre_T\tNSS_ave_pre_C\tNSS_ave_pre_G\tNSS_ave_firstBase_A\tNSS_ave_firstBase_T\tNSS_ave_firstBase_C\tNSS_ave_firstBase_G\n";
open SPECIES,$species_order_file or die "$0 fileOpenError: $species_order_file\n";
while(<SPECIES>){
	chomp;
	my $species=$_;
	if(exists $SS_feature_hash{$species} && exists $NSS_feature_hash{$species}){
		my $SS_out=$SS_feature_hash{$species};
		my $NSS_out=$NSS_feature_hash{$species};
		print FEATURE_OUT "$species\t$SS_out\t$NSS_out\n";
	}
	else{
		die "$species\n";
	}
}
close SPECIES;
close FEATURE_OUT;

##### FUNCTION
sub parse_basic{
	my($PmiREN2_species_file,$basic_file)=@_;
	my %spe_hash;
	my %spe_miRNA_hash;
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
		my($miRNA,$family,$mature,$species_id,$pre)=($line[13],$line[2],$line[18],$line[1],$line[9]);
		if(exists $spe_hash{$species_id}){
			$species=$spe_hash{$species_id};
		}
		else{die "$0 $species_id not in hash: spe_hash\n";}
		$spe_miRNA_hash{$species}{$miRNA}=$mature."\t".$pre;
	}
	close BASIC;
	return %spe_miRNA_hash;
}

sub base_composite{
	my($A,$T,$C,$G)=(0,0,0,0);
	my($seq)=@_;
	my @seq=split "",$seq;
	foreach my $base(@seq){
		if($base eq "A" || $base eq "a"){
			$A++;
		}
		elsif($base eq "T" || $base eq "t" || $base eq "U" || $base eq "u"){
			$T++;
		}
		elsif($base eq "C" || $base eq "c"){
			$C++;
		}
		elsif($base eq "G" || $base eq "g"){
			$G++;
		}
	}
	return $A,$T,$C,$G;
}

sub cal_feature_inSpecies{
	my %miRNA_num_hash; my %mature_len_hash; my %pre_len_hash; my %A_hash; my %T_hash; my %C_hash; my %G_hash; my %pre_A_hash; my %pre_T_hash; my %pre_C_hash; my %pre_G_hash; my %firstBase_hash;
	my %out_hash;
	my($infile)=@_;
	open INFILE,$infile or die "$0 fileOpenError: $infile\n";
	while(<INFILE>){
		chomp;
		my($species,$miRNA,$mature_len,$pre_len,$firstBase,$A,$T,$C,$G,$pre_A,$pre_T,$pre_C,$pre_G)=split "\t",$_;
		$miRNA_num_hash{$species}++;
		$mature_len_hash{$species}+=$mature_len;
		$pre_len_hash{$species}+=$pre_len;
		$A_hash{$species}+=$A;
		$T_hash{$species}+=$T;
		$C_hash{$species}+=$C;
		$G_hash{$species}+=$G;
		$pre_A_hash{$species}+=$pre_A;
		$pre_T_hash{$species}+=$pre_T;
		$pre_C_hash{$species}+=$pre_C;
		$pre_G_hash{$species}+=$pre_G;
		$firstBase_hash{$species}{$firstBase}++;
	}
	close INFILE;
	foreach my $species(sort keys %miRNA_num_hash){
		my $miRNA_num=$miRNA_num_hash{$species};
		if($miRNA_num != 1){
			my $ave_mature_len=$mature_len_hash{$species}/$miRNA_num;
			my $ave_pre_len=$pre_len_hash{$species}/$miRNA_num;
			my $ave_A=$A_hash{$species}/$miRNA_num;
			my $ave_T=$T_hash{$species}/$miRNA_num;
			my $ave_C=$C_hash{$species}/$miRNA_num;
			my $ave_G=$G_hash{$species}/$miRNA_num;
			my $ave_pre_A=$pre_A_hash{$species}/$miRNA_num;
			my $ave_pre_T=$pre_T_hash{$species}/$miRNA_num;
			my $ave_pre_C=$pre_C_hash{$species}/$miRNA_num;
			my $ave_pre_G=$pre_G_hash{$species}/$miRNA_num;
			my($ave_firstBase_A,$ave_firstBase_T,$ave_firstBase_C,$ave_firstBase_G)=(0,0,0,0);
			if(exists $firstBase_hash{$species}{"A"}){
				$ave_firstBase_A=$firstBase_hash{$species}{"A"}/$miRNA_num;
			}
			if(exists $firstBase_hash{$species}{"T"}){
				$ave_firstBase_T=$firstBase_hash{$species}{"T"}/$miRNA_num;
			}
			if(exists $firstBase_hash{$species}{"C"}){
				$ave_firstBase_C=$firstBase_hash{$species}{"C"}/$miRNA_num;
			}
			if(exists $firstBase_hash{$species}{"G"}){
				$ave_firstBase_G=$firstBase_hash{$species}{"G"}/$miRNA_num;
			}
			my $tmp1=$ave_mature_len."\t".$ave_pre_len;
			my $tmp2=$ave_A."\t".$ave_T."\t".$ave_C."\t".$ave_G;
			my $tmp3=$ave_pre_A."\t".$ave_pre_T."\t".$ave_pre_C."\t".$ave_pre_G;
			my $tmp4=$ave_firstBase_A."\t".$ave_firstBase_T."\t".$ave_firstBase_C."\t".$ave_firstBase_G;
			$out_hash{$species}=$tmp1."\t".$tmp2."\t".$tmp3."\t".$tmp4;
		}
		else{die "$0 miRNA_num in species: $species is zero\n";}
	}
	return %out_hash;
}







