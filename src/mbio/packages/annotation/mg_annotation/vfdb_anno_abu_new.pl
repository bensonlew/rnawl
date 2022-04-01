#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"c=s","pre=s","p=s","o=s","type=s","group=s");

my $usage = <<"USAGE";
        Program : $0
        Contact : shaohua.yuan\@majorbio.com
        Discription: Calculate vfdb abundance.
        Usage:perl $0 [options]
                -c     input gene_*_anno detail file
                -p     input gene profile file
                -o     output Dir
                -type  core or predict or all
                -group group file.

USAGE

die $usage  if (!($opts{c}&&$opts{o}&&$opts{p}&&$opts{type}));

`mkdir -p $opts{o}` unless -e ($opts{o});

my (%gene_profile,@sams,%profile_genes);
open INP,$opts{p} or die "can not open $opts{p}!\n";
while(<INP>){
    chomp;
    if (/^GeneID/){
    @sams = split /\t/;
    #### ɸ��head��GeneID�к�����total��
    shift @sams;
    #pop @sams;
    }else{
    my @tmp = split /\t/;
    my $gene = $tmp[0];
    $profile_genes{$gene} = 1;
    shift @tmp;
    my $tmp_profile = join("\t",@tmp);
    $gene_profile{$gene} = $tmp_profile;
    }
}
close INP;

my (%type_spe,%type_fun,%type_leve11,%type_level2,%type_gi,%type_vfs,%tatol_abu,%vf_l2,%datatype);
open INF,$opts{c} or die "can not open $opts{c}!\n";
while(<INF>){
    chomp;
    next if(/^#/);
    my @tmp = split /\t/;
    my $query = $tmp[0];
    if(exists $profile_genes{$query}){
        my $vfs = $tmp[4];
        my $gi = $tmp[1];
        my $l1 = $tmp[7];
        my $l2 = $tmp[8];
        my $spe = $tmp[6];
        my $fun = $tmp[5];
        if ($opts{type} =~ /all/){
            my $data_type = $tmp[-1];
            $type_spe{$data_type}{$vfs} = $spe;
            $type_fun{$data_type}{$vfs} = $fun;
            $type_leve11{$data_type}{$vfs} = $l1;
            $type_level2{$data_type}{$vfs} = $l2;
            $type_gi{$data_type}{$gi} .= ",".$query;
            $type_vfs{$data_type}{$vfs} .= ",".$query;
        }else{
            $type_spe{$vfs} = $spe;
            $type_fun{$vfs} = $fun;
            $type_leve11{$vfs} = $l1;
            $type_level2{$vfs} = $l2;
            $type_gi{$gi} .= ",".$query;
            $type_vfs{$vfs} .= ",".$query;
            }
        $vf_l2{$l1."_".$l2} .= ",".$query;
    }
}
close INF;


open OUP, "> $opts{o}/vfdb_level_pie.xls" or die "can not open $opts{o}/vfdb_level_pie.xls!\n";
my $samples = join("\t",@sams);

if (!($opts{type} =~ /all/)){
    open OUCG, "> $opts{o}/vfdb_$opts{type}_Gi_profile.xls" or die "can not open $opts{o}/vfdb_$opts{type}_Gi_profile.xls!\n";
    open OUCF, "> $opts{o}/vfdb_$opts{type}_VF_profile.xls" or die "can not open $opts{o}/vfdb_$opts{type}_VF_profile.xls!\n";
    print OUCG "#Gi_number\t","$samples\n";
    foreach my $eachgi (sort keys %type_gi){
        print OUCG "$eachgi\t";
        my $cgi_gene = $type_gi{$eachgi};
        my @genes = split(/,/,$cgi_gene);
        shift @genes;
        for(my $i = 0;$i < @sams ;$i++){
            my $cgi_abu;
            foreach my $gene (@genes){
                my @sams_profile = split(/\t/,$gene_profile{$gene});
                my $sam_profile = $sams_profile[$i];
                $cgi_abu += $sam_profile;
           }
        print OUCG $cgi_abu,"\t";
        }
        print OUCG "\n";
    }
    close OUCG;

    print OUCF "#VFs\t","$samples\tSpecies\tFunction\tLevel1\tLevel2\n";
    foreach my $vf (sort keys %type_vfs){
        print OUCF "$vf\t";
        my $cvf_gene = $type_vfs{$vf};
        my @genes = split(/,/,$cvf_gene);
        shift @genes;
        for(my $i = 0;$i < @sams ;$i++){
            my $cvf_abu;
            foreach my $gene (@genes){
                my @sams_profile = split(/\t/,$gene_profile{$gene});
                my $sam_profile = $sams_profile[$i];
                $cvf_abu += $sam_profile;
           }
        print OUCF $cvf_abu,"\t";
            }
        print OUCF "$type_spe{$vf}\t$type_fun{$vf}\t$type_leve11{$vf}\t$type_level2{$vf}\n";
    }
    close OUCF;
}else{
        print $opts{type},"\n";
        open OUAG, "> $opts{o}/vfdb_all_Gi_profile.xls" or die "can not open $opts{o}/vfdb_all_Gi_profile.xls!\n";
        open OUAF, "> $opts{o}/vfdb_all_VF_profile.xls" or die "can not open $opts{o}/vfdb_all_VF_profile.xls!\n";
        print OUAG "#Gi_number\tDatabase\t","$samples\n";
        print OUAF "#VFs\tDatabase\t","$samples\tSpecies\tFunction\tLevel1\tLevel2\n";
        foreach my $dtype(sort keys %type_gi){
        print $opts{type},"\n";
        open OUCG, "> $opts{o}/vfdb_".$dtype."_Gi_profile.xls" or die "can not open $opts{o}/vfdb_".$dtype."_Gi_profile.xls!\n";
        open OUCF, "> $opts{o}/vfdb_".$dtype."_VF_profile.xls" or die "can not open $opts{o}/vfdb_".$dtype."_VF_profile.xls!\n";
        print OUCG "#Gi_number\t","$samples\n";
        foreach my $eachgi (sort keys %{$type_gi{$dtype}}){
            print OUCG "$eachgi\t";
            print OUAG "$eachgi\t$dtype\t";
            my $cgi_gene = $type_gi{$dtype}{$eachgi};
            my @genes = split(/,/,$cgi_gene);
            shift @genes;
            for(my $i = 0;$i < @sams ;$i++){
                my $cgi_abu;
                foreach my $gene (@genes){
                    my @sams_profile = split(/\t/,$gene_profile{$gene});
                    my $sam_profile = $sams_profile[$i];
                    $cgi_abu += $sam_profile;
            }
            print OUCG "$cgi_abu\t";
            print OUAG "$cgi_abu\t";
            }
            print OUCG "\n";
            print OUAG "\n";
        }
        close OUCG;

        print OUCF "#VFs\t","$samples\tSpecies\tFunction\tLevel1\tLevel2\n";
        foreach my $vf (sort keys %{$type_vfs{$dtype}}){
            print OUCF "$vf\t";
            print OUAF "$vf\t$dtype\t";
            my $cvf_gene = $type_vfs{$dtype}{$vf};
            my @genes = split(/,/,$cvf_gene);
            shift @genes;
            for(my $i = 0;$i < @sams ;$i++){
                my $cvf_abu;
                foreach my $gene (@genes){
                    my @sams_profile = split(/\t/,$gene_profile{$gene});
                    my $sam_profile = $sams_profile[$i];
                    $cvf_abu += $sam_profile;
               }
            print OUCF $cvf_abu,"\t";
            print OUAF $cvf_abu,"\t";
                }
            print OUCF "$type_spe{$dtype}{$vf}\t$type_fun{$dtype}{$vf}\t$type_leve11{$dtype}{$vf}\t$type_level2{$dtype}{$vf}\n";
            print OUAF "$type_spe{$dtype}{$vf}\t$type_fun{$dtype}{$vf}\t$type_leve11{$dtype}{$vf}\t$type_level2{$dtype}{$vf}\n";
        }
        close OUCF;
      }
      close OUAG;
      close OUAF;
}

my (%group,%group_abu,%group_percet);
if ($opts{group}){
    open ING, $opts{group} or die "can not open $opts{group}!\n";
    while(<ING>){
        chomp;
        next if(/^#/);
        my @tmp = split /\t/;
        my $s_name = $tmp[0];
        my $g_name = $tmp[1];
        $group{$g_name} .= ",".$s_name;
    }
    close ING;
}

my %l2_abu;
if ($opts{group}){
    my $gs;
    foreach my $g(keys %group){
        $gs .= "\t".$g;
    }
    print OUP "Level1\tLevel2".$gs."\n";
}else{
    print OUP "Level1\tLevel2\t$samples\n";
}

#my %l2_abu;
my @l2_list;
my %l1_l2_hash;
foreach my $each_l2 (sort keys %vf_l2){
    my @l1_l2 = split(/_/,$each_l2);
    my $el1 = $l1_l2[0];
    my $el2 = $l1_l2[1];
    my $l2_gene = $vf_l2{$each_l2};
    my @genes = split(/,/,$l2_gene);
    shift @genes;
	my @el2_list = split(/;/,$el2);
	my @el1_list = split(/;/,$el1);
	for(my $j = 0;$j < @el2_list ;$j++){
		my $el2_s = $el2_list[$j];
		my $el1_s = $el1_list[$j];
		$el2_s =~ s/^\s+|\s+$//g;
		$el1_s =~ s/^\s+|\s+$//g;
		if(grep { $_ eq $el2_s } @l2_list){
			print "Have";
		}else{
			push @l2_list, $el2_s;
			$l1_l2_hash{$el2_s} = $el1_s;
		}
		for(my $i = 0;$i < @sams ;$i++){
			if(grep { $_ eq $el2_s } @l2_list){
				$l2_abu{$el2}{$sams[$i]} = $l2_abu{$el2}{$sams[$i]};
			}else{
				$l2_abu{$el2}{$sams[$i]} = 0;
			}
			foreach my $gene (@genes){
				my @sams_profile = split(/\t/,$gene_profile{$gene});
				my $sam_profile = $sams_profile[$i];
				$l2_abu{$el2_s}{$sams[$i]} += $sam_profile;
            }
		}

    #$tatol_abu{$el1}{$sams[$i]} += $l2_abu{$el2}{$sams[$i]};
    #print $el1,"\t",$l2_abu{$el2}{$sams[$i]},"\t",$tatol_abu{$el1}{$sams[$i]},"\n";
	}
    #$tatol_abu{$el1} += $l2_abu{$el2};
    #print $el1,"\t",$l2_abu{$el2},"\t",$tatol_abu{$el1},"\n";
}

foreach my $level2 (@l2_list){
	my $el2a = $level2;
	$el2a =~ s/^\s+|\s+$//g;
#	print $el2a, "\n";
    my $el1a = $l1_l2_hash{$el2a};
#	print $el1a, "\n";
	if($el1a ne "-"){
        print OUP "$el1a\t$el2a";
        if ($opts{group}){
            foreach my $g(keys %group){
                my @group_mems = split(/,/,$group{$g});
                for(my $i = 0;$i < @sams ;$i++){
                    if (grep /^$sams[$i]$/, @group_mems){
                        $group_abu{$g} += $l2_abu{$el2a}{$sams[$i]};
                        }
                }
                print OUP "\t$group_abu{$g}";
            }
        }else{
            for(my $m = 0;$m < @sams ;$m++){
                print OUP "\t$l2_abu{$el2a}{$sams[$m]}";
            }
        }
        print OUP "\n";
    }

}
close OUP;