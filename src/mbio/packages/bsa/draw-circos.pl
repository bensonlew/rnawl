#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($snp,$indel,$index,$index_column,$chrlist,$gff,$outdir,$windows,$gff_windows,$vcf);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(max min);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"windows:s"=>\$windows,
    "gwindows:s"=>\$gff_windows,
	"vcf:s"=>\$vcf,
    "index:s"=>\$index,
    "column:s"=>\$index_column,
	"gff:s"=>\$gff,
    "chrlist:s"=>\$chrlist,
	"outdir:s"=>\$outdir,
			) or &USAGE;
&USAGE unless ($windows and $vcf and $gff and $chrlist);
########
########work dir
########
$outdir||="./";
my $mkdir=1;
$mkdir=(mkdir "$outdir") if (!-d "$outdir");
$outdir=ABSOLUTE_DIR($outdir);
die "Error make dir $outdir" if ($mkdir == 0);
$mkdir=(mkdir "$outdir/draw.circos") if (!-d "$outdir/draw.circos");
die "Error make dir $outdir/draw.circos" if ($mkdir == 0);
$mkdir=(mkdir "$outdir/draw.circos/windows.file") if (!-d "$outdir/draw.circos/windows.file");
die "Error make dir $outdir/draw.circos/windows.file" if ($mkdir == 0);

open MC,">","$outdir/draw.circos/draw.conf";
open IDEO,">","$outdir/draw.circos/ideogram.conf";
open TIC,">","$outdir/draw.circos/ticks.conf";
open IM,">","$outdir/draw.circos/chromosomes.and.color.conf";
open HS,">","$outdir/draw.circos/housekeeping.conf";
my ($main_conf,$ideogram,$ticks,$chro,$housekeeping);
########
########chr and band file
########
open CHR,$chrlist;
open GFF,$gff;
open CHRBAND,">$outdir/draw.circos/chr.band.txt";
open CHRBSA,">$outdir/draw.circos/circos.chrlist";
print CHRBSA "[\n  ";
my %hash_chr_length;
my %hash_gff;
my %hash_chr_num;
my %hash_chr;
my $chr_type;
my @chrbsa;
while(<CHR>){
    #$_=~s/sca//g;
    $_=~s/[\n\r]//g;
	next if ($_ eq ""||/^$/||/^#/);
    my ($chr,$length)=split;
	if ($chr =~ /chr/) {
		$chr_type="chr";
	}else{
		$chr_type="sca";
	}
    $chr=~/(\D+)(\d+)/;
    my $chrnum=$2;
    $hash_chr{$chrnum}=1;
    $hash_chr_num{$chr}=1;
    $hash_chr_length{$chr}=$length;
#    print CHRBSA "  {\"id\":\"$chr\",\"label\":\"$chr\",\"color\":\"#996600\",\"len\":$length},\n";
	my $chrbsa="{\"id\":\"$chr\",\"label\":\"$chr\",\"color\":\"#996600\",\"len\":$length}";   
	push @chrbsa,"$chrbsa";
}
print CHRBSA join(",\n  ",@chrbsa);
print CHRBSA "]\n";
close CHRBSA;

my $chrnum=scalar keys %hash_chr_num;
$gff_windows||=$windows;
while(<GFF>){
    next if /^#/;
    $_=~s/[\n\r]//g;
    my (undef,undef,undef,undef,$chr,$start,@others)=split(/\t/,$_);
    my $win_num=int($start/$gff_windows)+1;
#    if ($chr =~ /chr/) {$chr=~s/chr//g;}
#	   if ($chr =~ /sca/) {$chr=~s/sca//g;}
	if(exists $hash_chr_num{$chr}){                ###########cuiqingmei
    $chr=~/(\D+)(\d+)/;
    $hash_gff{$1}{$2}{$win_num}++;   #chr= chr1/sca1;
	}
}
#foreach my $keys (sort {$hash_chr_num{$a}<=>$hash_chr_num{$b}} keys %hash_chr_num){
 #       print CHRBAND "chr"."\t"."-"."\t".$keys."\t".$hash_chr_num{$keys}."\t"."0"."\t"."$hash_chr_length{$keys}"."\t".$keys."\n";
#}
##draw colour of band
my %bands;
open CO,">$outdir/draw.circos/band.txt";
open COBSA,">$outdir/draw.circos/gene.num.csv";
print COBSA "Chr,start,end,GeneNum\n";
my @chr_windows;
my @draw_chr_cui;
foreach my $ctype (sort {$a cmp $b} keys %hash_gff){
	my $keys;
	my $chr_type=$ctype;
    	my $band=1;
	my $max=0;
	foreach $keys (sort {$a <=> $b} keys %{$hash_gff{$ctype}}){
	push @draw_chr_cui,$ctype.$keys;
    foreach my $num (sort keys %{$hash_gff{$ctype}{$keys}}){
        my $start=1+$gff_windows*($num-1);
        my $end=$gff_windows*$num; 
        my $snp=$hash_gff{$ctype}{$keys}{$num}/$gff_windows;
        push(@chr_windows,$snp);
        if($max < $num){
            $max=$num;
            }
        print CO "band\t$chr_type"."$keys\tband$band\tband$band\t$start\t$end\t$snp\n";
        print COBSA "$chr_type"."$keys,$start,$end,",$hash_gff{$ctype}{$keys}{$num},"\n"; ## gene.num.csv
        $band++;
    }
   }
    my $bandnum=$keys;
    $bandnum=~s/$chr_type//g; # cui
    print CHRBAND "chr"."\t"."-"."\t".$chr_type."".$keys."\t".$bandnum."\t"."0"."\t".$max*$gff_windows."\t".$chr_type."".$keys."\n";
}
close CO;
open BANDCOL,"<$outdir/draw.circos/band.txt";
my $band_windows=(max(@chr_windows)-min(@chr_windows))/10;
while(<BANDCOL>){
    $_=~s/[\n\r]//g;
    my @array=split;
    my $win_num=int($array[6]/$band_windows)+1;
    $array[5]=$array[5]-1;
    print CHRBAND "$array[0]\t$array[1]\t$array[2]\t$array[3]\t$array[4]\t$array[5]\tbandcol$win_num\n";
   # print CHRBAND join("\t",@array[1..($#array-1)])."\t"."bandcol".$win_num."\n";
}
#system("rm band.txt");
close CHR;
close GFF;
close CHRBAND;
########
########
########
#my $drawchr=join("",(1..$chrnum));
$chro.="
chromosomes_units = 1000000
chromosomes_display_default = yes
";
my @draw_chr=sort {$a <=> $b} keys %hash_chr;
# my $chr_draw=$chr_type.join(";$chr_type",@draw_chr);
my $chr_draw=join(";",@draw_chr_cui);
$chro.="
chromosomes=$chr_draw;
<colors>
";
my %chr_colour;
for(my $i=0;$i<@draw_chr;$i++){
	
    my $num=$i%6;
    if($i%6==0){
        $chro.="$chr_type$draw_chr[$i]* = col6\n";
        $chr_colour{$chr_type."".$draw_chr[$i]}="col6";
        }
    else{
        $chro.="$chr_type$draw_chr[$i]* = col$num\n";
        $chr_colour{$chr_type."".$draw_chr[$i]}="col".$num;
        }
    }
$chro.="
</colors>
";
#########
#########
#########
if (-e "$outdir/draw.circos/chr.band.txt"){
$main_conf.="
karyotype = $outdir/draw.circos/chr.band.txt

<<include $outdir/draw.circos/ideogram.conf>>
<<include $outdir/draw.circos/ticks.conf>>
<<include $outdir/draw.circos/chromosomes.and.color.conf>>
<image>
<<include etc/image.conf>>
radius* = 2000
</image>
";
}
########
########
########
$ideogram="
<ideogram>
<spacing>
default = 0.005r
</spacing>
radius    = 0.9r
thickness = 100p
fill      = yes
fill_color = black
#stroke_thickness = 2
#stroke_color     = black
show_label       = yes
label_font       = default 
label_radius     = 1r + 75p
label_size       = 30
label_parallel   = yes
label_case     = upper 
label_format   = eval(sprintf(\"$chr_type%s\",var(label)))
show_bands            = yes
fill_bands            = yes
band_stroke_thickness = 2
band_stroke_color     = white
band_transparency     = 1
</ideogram>
";
#########
#########
#########
$ticks="
show_ticks          = no
show_tick_labels    = no
<ticks>
radius           = 1r
color            = black
thickness        = 2p
multiplier       = 1e-6
format           = %d
<tick>
spacing        = 5u
size           = 10p
</tick>
<tick>
spacing        = 25u
size           = 15p
show_label     = yes
label_size     = 20p
label_offset   = 10p
format         = %d
</tick>

</ticks>
";
$housekeeping.="
anglestep       = 0.5
minslicestep    = 10
beziersamples   = 40
debug           = no
warnings        = no
imagemap        = no
paranoid        = no
units_ok        = bupr
units_nounit    = n
file_delim = \\s
file_delim_collapse = yes
list_record_delim = \\s*[;,]\\s*
list_field_delim  = \\s*[:=]\\s*
options_record_delim = [,;]
options_field_delim  = =
skip_missing_expression_vars = no
legacy_underline_expression_syntax = no
svg_font_scale = 1.3
sup_baseline_shift = 40
sub_baseline_shift = -40
sup_fontsize = 90
sub_fontsize = 90
default_font   = default
default_font_name  = Arial
default_font_color = black
default_color  = black
<guides>
thickness      = 1
size           = 5
type           = outline
<object>
all            = no
ideogram       = no
ideogram_label = no
</object>
<color>
default = lblue
text    = red
</color>
</guides>
debug_group = summary,output
debug_auto_timer_report = 30
debug_word_separator = \" \"
debug_undef_text     = _undef_
debug_empty_text     = _emptylist_
debug_validate       = yes
debug_output_tidy    = no
text_pixel_subsampling = 1
text_snuggle_method    = array
restrict_parameter_names = no
case_sensitive_parameter_names = no
calculate_track_statistics = yes
color_cache_static = yes
color_cache_file   = circos.colorlist
color_lists_use    = yes
memoize = yes
quit_on_dump = yes
offsets = 0,0
max_ticks            = 5000
max_ideograms        = 200
max_links            = 25000000000000000000
max_points_per_track = 25000000000000000000
undefined_ideogram = skip
relative_scale_iterations = 10
relative_scale_spacing    = mode
data_out_of_range = trim,warn
track_defaults = etc/tracks
round_brush_use           = yes
round_brush_min_thickness = 5
anti_aliasing = yes
housekeeping = yes
auto_eval = no
";
#########
$main_conf.="
<plots>
<backgrounds>
show  = data
<background>
color = vvlgrey
y0    = 1.0r
y1    = 0r
</background>
</backgrounds>
";
#    open IN,"<$vcf";
########################################
if ($vcf =~ /gz$/) {
        open IN,"gunzip -dc $vcf|";
}else{
        open IN,$vcf;
}
########################################
    open SNP,">$outdir/draw.circos/windows.file/snp.win.txt";
    open SNPBSA,">$outdir/draw.circos/snp.win.csv";
print SNPBSA "Chr,start,end,SNPNum\n";
	open INDEL,">$outdir/draw.circos/windows.file/indel.win.txt";
	open INDELBSA,">$outdir/draw.circos/indel.win.csv";
print INDELBSA "Chr,start,end,INDELNum\n";
	my %max;
	slide($vcf,\%hash_chr_num,\%max);
close IN;
    $main_conf.="
		######plot snp
		<plot>
		type = line
		max_gap = 1u
		file    = $outdir/draw.circos/windows.file/snp.win.txt
		color   = white
		min     = 0
		max     = $max{snp} #0.015
		r0      = 0.85r#1.075r
		r1      = 0.95r#1.125r
		thickness = 1
		#fill_color = vdyellow
		</plot>
		######
		";
		$main_conf.="
		<plot>
			type=scatter
			file= $outdir/draw.circos/windows.file/indel.win.txt
			#fill_color=green
			stroke_color=blue
			glyph=rectangle#circle
			glyph_size=10
			max=$max{indel}#0.013
			min=0
			r1=0.75r
			r0=0.65r
		</plot>
		";

if($index){######plot sv
    my $file_name=basename($index);
    open IN,"<$index";
    open OUT,">$outdir/draw.circos/windows.file/$file_name.win.txt";
    open OUTBSA,">$outdir/draw.circos/sliding.win.csv";
    print OUTBSA "Chr,start,end,index\n";
    my (@max)=slide_index($index,\%hash_chr_num);
    $main_conf.="
<plot>
type = line
max_gap = 1u
file = $outdir/draw.circos/windows.file/$file_name.win.txt
color = white
min = $max[1]
max = $max[0]#0.0001
r0         = 0.45r
r1         = 0.55r
thickness = 0.1
#fill_under = yes
#fill_color = blue
#orientation      = out
</plot>
";
}

$main_conf.="
</plots>
<colors>
<<include $outdir/draw.circos/colors.conf>>
</colors>
<<include etc/colors_fonts_patterns.conf>>
<<include $outdir/draw.circos/housekeeping.conf>>
";

print MC $main_conf;
print IDEO $ideogram;
print TIC $ticks;
print IM $chro;
print HS $housekeeping;

close MC;
close IDEO;
close TIC;
close IM;
open Color,">$outdir/draw.circos/colors.conf";
print Color "col1=rgb(255,0,0)\n";
print Color "col2=rgb(0,255,0)\n";
print Color "col3=rgb(0,0,255)\n";
print Color "col4=rgb(0,255,255)\n";
print Color "col5=rgb(255,0,255)\n";
print Color "col6=rgb(255,165,0)\n";
print Color "white=rgb(0,0,0)\n";
print Color "black=rgb(205,205,205)\n";
print Color "bandcol1=rgb(195,195,195)\n";
print Color "bandcol2=rgb(185,185,185)\n";
print Color "bandcol3=rgb(175,175,175)\n";
print Color "bandcol4=rgb(165,165,165)\n";
print Color "bandcol5=rgb(155,155,155)\n";
print Color "bandcol6=rgb(145,145,145)\n";
print Color "bandcol7=rgb(135,135,135)\n";
print Color "bandcol8=rgb(125,125,125)\n";
print Color "bandcol9=rgb(115,115,115)\n";
print Color "bandcol10=rgb(105,105,105)\n";
close Color;
system("circos -conf $outdir/draw.circos/draw.conf -outputfile circos -outputdir $outdir/ ");

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub slide{
    my ($file,$chr,$max)=@_;
    my %hash;
    my $win=$windows;
    my %hash_max;
	$$max{snp}||=0;
	$$max{indel}||=0;
    while(<IN>){
        #$_=~s/sca/chr/g;
        $_=~s/[\n\r]//g;
        my @array=split;

        next if /^#/;
        my $win_num=int($array[1]/$win)+1;
		next if (!exists $$chr{$array[0]});
		my $lenr=length($array[3]);
		my @all=split(/,/,$array[4]);
		my $lena=0;
		for (my $i=0;$i<@all;$i++) {
			if (length($all[$i]) > $lena) {
				$lena=length($all[$i]);
			}
		}
		if ($lena == $lenr) {
			$hash{$array[0]}{$win_num}{snp}++;
		}else{
			$hash{$array[0]}{$win_num}{indel}++;
		}
    }
    foreach my $keys (sort keys %hash){
   # print "$keys\n";
        foreach my $num (sort keys %{$hash{$keys}}){
			$hash{$keys}{$num}{snp}||=0;
			$hash{$keys}{$num}{indel}||=0;
            my $chr=$keys;
            my $start=1+$win*($num-1);
            my $end=$win*$num;
            my $snp=$hash{$keys}{$num}{snp}/$win;
			my $indel=$hash{$keys}{$num}{indel}/$win;
			if ($snp > $$max{snp}) {
				$$max{snp}=$snp;
			}
			if ($indel > $$max{indel}) {
				$$max{indel}=$indel;
			}
			print SNP "$chr\t$start\t$end\t$snp\tfill_color=$chr_colour{$chr}\n";
			print SNPBSA "$chr,$start,$end,",$hash{$keys}{$num}{snp},"\n";
			print INDEL "$chr\t$start\t$end\t$indel\tfill_color=$chr_colour{$chr}\n";
			print INDELBSA "$chr,$start,$end,",$hash{$keys}{$num}{indel},"\n";
        
    }
    }
}
sub slide_index{
    my ($file,$chr)=@_;
    my @max;
    my %hash;
    my $line=0;
    my $count=0;
    my $chra="";
    while(<IN>){
        $_=~s/[\n\r]//g;
        if($count >=1 ){
            my @array=split;
            my $deltaSNP=$array[$index_column-1];
#            $array[0]=~/\"(\D+)(\d+)\"/;
#            my $chr=$chr_type."".$2;a
            my $chr;
            $array[0]=~s/"//g;
	    if(exists $hash_chr_num{$array[0]}){
	            $chr=$array[0];
	    }
            $hash{$deltaSNP}=1;
            if($array[1]=~/(\d+)e\+(\d+)/){
                $array[1]=$1*(10**$2);
                }
            if($array[2]=~/(\d+)e\+(\d+)/){
                $array[2]=$1*(10**$2);
            }
            next if (!defined $chr);
            if($chr eq $chra){
                if($line==0){
                    print OUT "$chr\t$array[1]\t$array[2]\t$deltaSNP\tfill_color=$chr_colour{$chr}\n";
                    print OUTBSA "$chr,$array[1],$array[2],$deltaSNP\n";
                }else{
                    my $newstart=$array[2]-10000;
                    print OUT "$chr\t$newstart\t$array[2]\t$deltaSNP\tfill_color=$chr_colour{$chr}\n";
                    print OUTBSA "$chr,$newstart,$array[2],$deltaSNP\n";
                    }
                    $line++;
            }else{
                $line=0;
                if($line==0){
                    print OUT "$chr\t$array[1]\t$array[2]\t$deltaSNP\tfill_color=$chr_colour{$chr}\n";
                    print OUTBSA "$chr,$array[1],$array[2],$deltaSNP\n";
                }else{
                    my $newstart=$array[2]-10000;
                    print OUT "$chr\t$newstart\t$array[2]\t$deltaSNP\tfill_color=$chr_colour{$chr}\n";
                    print OUTBSA "$chr,$newstart,$array[2],$deltaSNP\n";
                    }
                    $chra=$chr;
                    $line++;
                }
        }
        $count++;
    }
    $max[0]=(sort{$b<=>$a} keys %hash)[0];
    $max[1]=(sort{$a<=>$b} keys %hash)[0];
    return @max;
}
sub ABSOLUTE_DIR 
{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir \n$in";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub USAGE {#
        my $usage=<<"USAGE";
Contact:        nanshan.yang\@majorbio.com;
Script:			$Script
Description:
	draw circos
	eg:
	perl $Script 

Usage:
  Options:
  --windows windows numbers
  --gwindows gff file windows numbers
  --snp	<file> slide windows of snp
  --indel <file> slide windows of indel
  --index   <file> mutmap index file 
  --column  <number> the column of the index file
  --chrlist input chrome number
  --gff gff file
  -h         Help

USAGE
        print $usage;
        exit;
}
