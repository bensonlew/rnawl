#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($snp,$indel,$sv,$cnv,$chrlist,$gff,$outdir,$windows);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(max min);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"windows:s"=>\$windows,
	"snp:s"=>\$snp,
    "indel:s"=>\$indel,
    "sv:s" =>\$sv,
    "cnv:s"=>\$cnv,
	"gff:s"=>\$gff,
    "chrlist:s"=>\$chrlist,
	"outdir:s"=>\$outdir,
			) or &USAGE;
&USAGE unless ($windows and $snp and $gff and $chrlist);
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
my %hash_chr_length;
my %hash_gff;
my %hash_chr_num;
my $n=0;
my %hash_chr;
while(<CHR>){
    $_=~s/[\n\r]//g;
	next if ($_ eq ""||/^$/||/^#/);
    my ($chr,$length)=split;
    $n++;
    $hash_chr_num{$chr}=$n;
    $hash_chr{$n}=$chr;
    $hash_chr_length{$chr}=$length;
}
my $chrnum=scalar keys %hash_chr_length;
my $gff_windows=$windows;

while(<GFF>){
    next if /^#/;
    $_=~s/[\n\r]//g;
    my ($chr,undef,$type,$start,@others)=split(/\t/,$_);
    my $win_num=int($start/$gff_windows)+1;
    if (exists $hash_chr_num{$chr} and ($type=~/gene/ || ($type !~ /CDS/ && $type !~/exon/)) ){        
        $hash_gff{$hash_chr_num{$chr}}{$win_num}++;
        }
}

my %bands;
open CO,">$outdir/draw.circos/band.txt";
my @chr_windows;
foreach my $keys (sort {$hash_chr_num{$a}<=>$hash_chr_num{$b}}keys %hash_chr_num){
    my $band=1;
    my $max=0;
    my $chr=$keys;
	if ($hash_chr_length{$chr} < $gff_windows) {
		$hash_chr_length{$chr}=$gff_windows;
	}
	if (!exists $hash_gff{$hash_chr_num{$keys}}) {
		my $bandnum=$hash_chr_num{$keys};
		print CHRBAND "chr"."\t"."-"."\t".$chr."\t".$chr."\t"."0"."\t".int($hash_chr_length{$chr}/$gff_windows+1)*$gff_windows."\t"."".$keys."\n";
		next;
	}
    foreach my $num (sort keys %{$hash_gff{$hash_chr_num{$keys}}}){
        my $start=1+$gff_windows*($num-1);
        my $end=$gff_windows*$num; 
        my $snp=$hash_gff{$hash_chr_num{$keys}}{$num}/$gff_windows;
        push(@chr_windows,$snp);
        if($max < $num){
            $max=$num;
            }
        print CO "band\t$chr\tband$band\tband$band\t$start\t$end\t$snp\n";
        $band++;
    }
    my $bandnum=$hash_chr_num{$keys};
    print CHRBAND "chr"."\t"."-"."\t".$chr."\t".$chr."\t"."0"."\t".int($hash_chr_length{$chr}/$gff_windows+1)*$gff_windows."\t"."".$keys."\n";
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
}

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
my @draw_chr=sort {$hash_chr_num{$a} <=> $hash_chr_num{$b}} keys %hash_chr_num;
my $chr_draw=join(";",@draw_chr);
$chro.="
chromosomes=$chr_draw;
<colors>
";
my %chr_colour;
for(my $i=0;$i<@draw_chr;$i++){
	
    my $num=$i%6;
    if($i%6==0){
        $chro.="$draw_chr[$i]* = col6\n";
        $chr_colour{$draw_chr[$i]}="col6";
        }
    else{
        $chro.="$draw_chr[$i]* = col$num\n";
        $chr_colour{$draw_chr[$i]}="col".$num;
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
#label_format   = eval(sprintf(\"%s\",var(label)))
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
print $snp;
if ($snp){######plot snp
    #print $snp;
    my $file_name=(split/\//,$snp)[-1];
	open IN,$snp;
	if ($snp=~/.gz$/) {
		close IN;
		open IN,"gunzip -c $snp|";
	}
    open OUT,">$outdir/draw.circos/windows.file/$file_name.win.txt";
    my $max=slide($snp,\%hash_chr_num);
    $main_conf.="
######plot snp
<plot>
type = line
max_gap = 1u
file    = $outdir/draw.circos/windows.file/$file_name.win.txt
color   = white
min     = 0
max     = $max#0.015
r0      = 0.85r#1.075r
r1      = 0.95r#1.125r
thickness = 1
#fill_color = vdyellow
</plot>
######
";
}
print $indel,"\n";
if($indel){######plot indel
    my $file_name=(split/\//,$indel)[-1];
	open IN,$indel;
	if ($indel=~/.gz$/) {
		close IN;
		open IN,"gunzip -c $indel|";
	}
    open OUT,">$outdir/draw.circos/windows.file/$file_name.win.txt";
    my $max=slide($indel,\%hash_chr_num);
    $main_conf.="
<plot>
type=scatter
file= $outdir/draw.circos/windows.file/$file_name.win.txt
#fill_color=green
stroke_color=blue
glyph=rectangle#circle
glyph_size=10
max=$max#0.013
min=0
r1=0.75r
r0=0.65r
</plot>
";
}
if($snp and $indel and $sv){######plot sv
    my $file_name=(split/\//,$sv)[-1];
    open IN,"<$sv";
    open OUT,">$outdir/draw.circos/windows.file/$file_name.win.txt";
    my $max=slide($sv,\%hash_chr_num);
    $main_conf.="
<plot>
type = histogram
file = $outdir/draw.circos/windows.file/$file_name.win.txt
color = white
min = 0
max = $max#0.0001
r0         = 0.45r
r1         = 0.55r
thickness = 2
fill_under = yes
#fill_color = blue
orientation      = out
</plot>
";
}
if($snp and $indel and $sv and $cnv){######plot cnv
    my $file_name=(split/\//,$cnv)[-1];
    open IN,"<$cnv";
    open OUT,">$outdir/draw.circos/windows.file/$file_name.win.txt";
    my $max= slide($cnv,\%hash_chr_num);
    $main_conf.="
<plot>
type = scatter
file = $outdir/draw.circos/windows.file/$file_name.win.txt
#fill_color=blue
stroke_color=black
glyph=circle
glyph_size=10
min     = 0
max     = $max#0.001
r0   = 0.25r
r1   = 0.35r
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
    my ($file,$chr)=@_;
    my $max;
    my %hash;
    my $win=$windows;
    my %hash_max;
    while(<IN>){
        #$_=~s/sca/chr/g;
        $_=~s/[\n\r]//g;
        my @array=split;
        next if /^#/;
        my $win_num=int($array[1]/$win)+1;
		next if (!exists $$chr{$array[0]});
        $hash{$array[0]}{$win_num}++;
    }
    foreach my $keys (sort keys %hash){
   # print "$keys\n";
        foreach my $num (sort keys %{$hash{$keys}}){
            my $chr=$keys;
            my $start=1+$win*($num-1);
            my $end=$win*$num;
            my $snp=$hash{$keys}{$num}/$win;
			$hash_max{$snp}=1;
            print OUT "$chr\t$start\t$end\t$snp\tfill_color=$chr_colour{$chr}\n";
        
    }
    }
    $max= (sort{$b<=>$a} keys %hash_max)[0];
    return $max;
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
  --snp	<file> slide windows of snp
  --indel <file> slide windows of indel
  --sv <file> slide windows of sv
  --cnv <file> slide windows of cnv
  --chrlist input chrome number
  --gff gff file
  --outdir
  -h         Help

USAGE
        print $usage;
        exit;
}
