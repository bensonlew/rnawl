#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/3/28 16:44
@file    : circos_config.py
"""

ideogram = '''
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
'''

colors = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
          'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chrX', 'chr24',
          'chrY', 'chrM', 'chr0', 'chrUn', 'chrNA']

colors_string = '''
# chromosome color map (UCSC) 

chr1 = 153,102,0
chr2 = 102,102,0
chr3 = 153,153,30
chr4 = 204,0,0
chr5 = 255,0,0
chr6 = 255,0,204
chr7 = 255,204,204
chr8 = 255,153,0
chr9 = 255,204,0
chr10 = 255,255,0
chr11 = 204,255,0
chr12 = 0,255,0
chr13 = 53,128,0
chr14 = 0,0,204
chr15 = 102,153,255
chr16 = 153,204,255
chr17 = 0,255,255
chr18 = 204,255,255
chr19 = 153,0,204
chr20 = 204,51,255
chr21 = 204,153,255
chr22 = 102,102,102
chr23 = 153,153,153
chrX  = 153,153,153
chr24 = 204,204,204
chrY = 204,204,204
chrM = 204,204,153
chr0 = 204,204,153
chrUn = 121,204,61
chrNA = 255,255,255

white=rgb(0,0,0)
black=rgb(205,205,205)
bandcol1=rgb(195,195,195)
bandcol2=rgb(185,185,185)
bandcol3=rgb(175,175,175)
bandcol4=rgb(165,165,165)
bandcol5=rgb(155,155,155)
bandcol6=rgb(145,145,145)
bandcol7=rgb(135,135,135)
bandcol8=rgb(125,125,125)
bandcol9=rgb(115,115,115)
bandcol10=rgb(105,105,105)

item_grep=128,128,128
item_red=255,0,0
item_blue=0,0,255
item_chen=255,102,51
item_green=0,255,0

'''

chromosomes_and_color = '''
chromosomes_units = 1000000
chromosomes_display_default = yes

chromosomes = {chromosomes}

<colors>
{colors}
</colors>
'''

housekeeping = '''
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

debug_group                 = summary,output
debug_auto_timer_report     = 30
debug_word_separator        = \" \"
debug_undef_text            = _undef_
debug_empty_text            = _emptylist_
debug_validate              = yes
debug_output_tidy           = no
text_pixel_subsampling      = 1
text_snuggle_method         = array
restrict_parameter_names    = no
case_sensitive_parameter_names = no
calculate_track_statistics  = yes
color_cache_static          = yes
color_cache_file            = circos.colorlist
color_lists_use             = yes
memoize                     = yes
quit_on_dump                = yes
offsets                     = 0,0
max_ticks                   = 5000
max_ideograms               = 200
max_links                   = 25000000000000000000
max_points_per_track        = 25000000000000000000
undefined_ideogram          = skip
relative_scale_iterations   = 10
relative_scale_spacing      = mode
data_out_of_range           = trim,warn
track_defaults              = etc/tracks
round_brush_use             = yes
round_brush_min_thickness   = 5
anti_aliasing               = yes
housekeeping                = yes
auto_eval                   = no
'''

ticks = '''
show_ticks       = no
show_tick_labels = no
<ticks>
radius           = 1r
color            = black
thickness        = 2p
multiplier       = 1e-6
format           = %d
<tick>
spacing          = 5u
size             = 10p
</tick>
<tick>
spacing          = 25u
size             = 15p
show_label       = yes
label_size       = 20p
label_offset     = 10p
format           = %d
</tick>
</ticks>
'''

histogram = '''
<plot>
type        = histogram
file        = {file}
color       = white
min         = {min}
max         = {max}
r0          = {r0}
r1          = {r1}
thickness   = {thickness}
fill_under  = yes
# fill_color = {fill_color}
fill        = yes
extend_bin = no
# orientation = {orientation}
</plot>
'''

links = '''
<link>
file = {file}
radius = {radius}
bezier_radius = 0r
color = {color}
thickness = {thickness}
</link>
'''

main_conf = '''
karyotype = {karyotype}

<<include {ideogram_conf}>>
<<include {ticks_conf}>>
<<include {chromosomes_and_color_conf}>>

<image>
<<include etc/image.conf>>
radius* = 2000
</image>

<plots>

<backgrounds>
show  = data
<background>
color = vvlgrey
y0    = 1.0r
y1    = 0r
</background>
</backgrounds>

# plot 
{plots}

</plots>

<links>

{links}

</links>

<colors>
<<include {colors_conf}>>
</colors>

<<include etc/colors_fonts_patterns.conf>>
<<include {housekeeping_conf}>>
'''

chr_colors = [
    '153,102,0', '102,102,0', '153,153,30', '204,0,0', '255,0,0', '255,0,204', '255,204,204', '255,153,0',
    '255,204,0', '255,255,0', '204,255,0', '0,255,0', '53,128,0', '0,0,204', '102,153,255', '153,204,255',
    '0,255,255', '204,255,255', '153,0,204', '204,51,255', '204,153,255', '204,204,153', '121,204,61'
]

band_colors = ['bandcol1', 'bandcol2', 'bandcol3', 'bandcol4', 'bandcol5',
               'bandcol6', 'bandcol7', 'bandcol8', 'bandcol9', 'bandcol10']
chr_data = {
    'chr': 'chr\t-\t{chrom}\t{chrom_lable}\t0\t{end:.0f}\t{color}\n',
    'gene': 'band\t{chrom}\t{band}\t{band}\t{start:.0f}\t{end:.0f}\t{color}\n'
}

mrna_histo_demo = '{chrom}\t{start:.0f}\t{end:.0f}\t{value}\tfill_color={color}\n'
# segdup00011 hs1 71096 76975 thickness=5 segdup00011 hs1 388076 393885 thickness=5
cir_line_demo = '{chrom}\t{start:.0f}\t{end:.0f}\t' \
                '{t_chrom}\t{t_start:.0f}\t{t_end:.0f}\tcolor={color}\n'
