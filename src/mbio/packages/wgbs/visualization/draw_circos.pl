#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($pos, $neg, $label, $chrlist, $chr_band, $gff,$outdir,$windows);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(max min);
my $version="1.0.0";
GetOptions(
    "help|?" =>\&USAGE,
    "windows:s"=>\$windows,
    "pos:s"=>\$pos,
    "neg:s"=>\$neg,
    "label:s"=>\$label,
    "gff:s"=>\$gff,
    "chr_band:s"=>\$chr_band,
    "chrlist:s"=>\$chrlist,
    "outdir:s"=>\$outdir,
			) or &USAGE;
&USAGE unless ($pos and $neg  and $chrlist);
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
open BREAK,">","$outdir/draw.circos/break.conf";
open TIC,">","$outdir/draw.circos/ticks.conf";
open IM,">","$outdir/draw.circos/chromosomes.and.color.conf";
open HS,">","$outdir/draw.circos/housekeeping.conf";
my ($main_conf,$ideogram,$ticks,$chro,$housekeeping, $break);
########
########chr and band file
########
open CHR,$chrlist;

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


close CHR;


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
        $chro.="$draw_chr[$i]* = white\n";
        $chr_colour{$draw_chr[$i]}="white";
        }
    else{
        $chro.="$draw_chr[$i]* = white\n";
        $chr_colour{$draw_chr[$i]}="white";
        }
    }
$chro.="
</colors>
";
if ($chr_band){
    system("cp $chr_band $outdir/draw.circos/chr.band.txt");
}
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
$break="
axis_brea
axis_break         = yes
axis_break_style   = 2

<break_style 1>
stroke_color     = black
fill_color       = blue
thickness        = 0.25r
stroke_thickness = 2
</break_style>

<break_style 2>
stroke_color     = black
stroke_thickness = 2
thickness        = 2r
</break_style>
";

########
########
### <<include $outdir/draw.circos/break.conf>>
########
$ideogram="
<ideogram>
<spacing>
default = 0.005r

</spacing>
radius    = 0.6r
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
label_format   = eval(sprintf(\"%s\",var(label)))
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
show_ticks          = yes
show_tick_labels    = yes
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
color = white
y0    = 1.0r
y1    = 0r
</background>
</backgrounds>
";

=pod
<axes>
<axis>
axis = yes
color     = lgrey
thickness = 1
spacing   = 0.3r
show_label     = yes
position_skip = 0.5
show_ticks       = yes
show_tick_labels = yes
</axis>
</axes>
=cut
    
print $pos;
if ($pos){
    ######plot snp
    #print $snp;
    system("cp $pos $outdir/draw.circos/windows.file/pos.win.txt");
    $main_conf.="
######plot pos histogram
<plot>
type = histogram
file    = $outdir/draw.circos/windows.file/pos.win.txt
color   = red
r0      = 0.6r
r1      = 0.85r
thickness = 2
fill_under = yes
orientation  = out
fill_color = red
min = 0
max = 0.4

</plot>
######
";
}

print $neg;
if ($neg){
    ######plot snp
    #print $snp;
    system("cp $neg $outdir/draw.circos/windows.file/neg.win.txt");
    $main_conf.="
######plot neg histogram
<plot>
type = histogram
file    = $outdir/draw.circos/windows.file/neg.win.txt
color   = blue
r0      = 0.35r
r1      = 0.6r
thickness = 2
fill_under = yes
orientation  = in
fill_color = blue
min = 0
max = 0.4

</plot>
######
";
}

print $label;
if ($label){
    ######plot label
    #print $label;
    system("cp $label $outdir/draw.circos/windows.file/label.txt");
    $main_conf.="
######plot label histogram
<plot>
type             = text
color            = black
file             = $outdir/draw.circos/windows.file/label.txt
r0 = 1.0r
r1 = 1.2r
show_links     = yes
link_dims      = 4p,4p,8p,4p,4p
link_thickness = 2p
link_color     = red
label_size   = 24p
fill_under = yes
label_font   = condensed
padding  = 0p
rpadding = 0p
</plot>
######
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
print BREAK $break;
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
  --pos <file> pos mapping windows
  --neg <file> neg mapping windows
  --chrlist input chrome number
  --gff gff file
  --outdir
  -h         Help

USAGE
        print $usage;
        exit;
}
