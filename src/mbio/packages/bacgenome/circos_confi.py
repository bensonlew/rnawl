# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'

import os
import argparse

def cal_space(gsize):
    spc1  = gsize/20
    if spc1 < 1000  :
        spc = 1000
        spc2 = 200
    elif 1000 < spc1 < 3001 :
        spc = 1000
        spc2 = 200
    elif 3000 < spc1 < 20000 :
        spc = 10000
        spc2 = 2000
    else:
        spc = 100000
        spc2 = 20000

    return spc,spc2




def get_genome_size(infile):  #get from karyotype.txt : chr - p 1 1 194272 deepskyblue
    with open(infile) as fr:
        line = fr.readline()
        try :
            gsize = int(line.split()[-2])
            return gsize
        except:
            return None

def write_txt(output, kary, cog, temp, an_cog, p_gc_c, n_gc_c, p_gc_s, n_gc_s, file=None, type=1, circle=None):
    if circle == "Circular":
        if file:
            print(file)
        with open(output, 'w') as o:
            if int(type) == 1:
                multiplier = "1e-6"
                spacing1 = "20000b"
                spacing2 = "100000b"
            else:
                multiplier = "1e-3"
                spacing1 = "200b"
                spacing2 = "1000b"
                gsize = get_genome_size(kary)
                if gsize:
                    big, small = cal_space(gsize)
                    spacing2 = str(big) + 'b'
                    spacing1 = str(small) + 'b'

            o.write('<<include etc/colors_fonts_patterns.conf>>\n<image>\n<<include etc/image.conf>>\n</image>\n\n')
            o.write('karyotype = ' + kary + '\n')
            o.write(
                'chromosomes_units = 100000\nchromosomes_display_default = yes\n\n<ideogram>\n<spacing>\n'
                'default = 0.005r\n</spacing>\nradius = 0.80r\nthickness = 1p\nfill = yes\nfill_color = deepskyblue\nstroke_color = black\n'  # thickness = 6p 改成 1p
                'stroke_thickness = 1p\nshow_label = no\nlabel_font = light\nlabel_radius = 1r + 110p\nlabel_size = 30\nlabel_parallel = no\n'
                '</ideogram>\n\nshow_ticks = yes\nshow_tick_labels = yes\n<ticks>\nskip_first_label = no\nskip_last_label = no\n'
                'radius = dims(ideogram,radius_outer)\ncolor = black\nthickness = 2p\nsize = 30p\nmultiplier = ' + multiplier + '\nformat = %.2f'
                                                                                                                                '\n<tick>\nspacing = ' + spacing1 + '\nsize = 10p\nshow_label = no\nthickness = 3p\n</tick>\n<tick>\nspacing = ' + spacing2 + '\nsize = 20p\nshow_label = yes\n'
                                                                                                                                                                                                                                                              'label_size = 12p\nlabel_offset = 10p\nformat = %.2f\n</tick>\n</ticks>\n\n<highlights>\nz = 0\n<highlight>\n')
            # label_size 25 改成12
            o.write('file = ' + cog + '\n')
            o.write('r0 = 0.89r\nr1 = 0.99r\n</highlight>\n<highlight>\n')
            # o.write('file = ' + temp + '\n')
            # o.write('r0 = 0.89r\nr1 = 0.8901r\nfill_color = black\n</highlight>\n<highlight>\n')  #r1 = 0.891r  r1 = 0.8901r
            # o.write('r0 = 0.89r\nr1 = 0.891r\nstroke_color = black\nstroke_thickness = 1p\n</highlight>\n<highlight>\n')
            o.write('file = ' + an_cog + '\n')
            o.write('r0 = 0.79r\nr1 = 0.89r\n</highlight>\n')
            if file:
                file_list = file.split(',')
                div = (0.78 - 0.60) / len(file_list)
                n = 0
                for i in file_list:
                    o.write('<highlight>\nfile = ' + i + '\nr0 = ' + str(0.78 - n * div) + 'r\n' + 'r1 = ' + str(
                        0.78 - (n + 1) * div) + 'r\n</highlight>\n')
                    n = n + 1
            o.write('</highlights>\n\n')
            o.write('<plots>\ntype = line\nthickness = 1p\n<plot>\nz = 2\nmax_gap = 0u\n')
            o.write('file = ' + p_gc_c + '\n')
            o.write(
                'color = red\nfill_color = red\nr0 = 0.48r\nr1 = 0.58r\norientation = out\n</plot>\n<plot>\nz = 2\nmax_gap = 0u\n')
            o.write('file = ' + n_gc_c + '\n')
            o.write(
                'color = blue\nfill_color = blue\nr0 = 0.38r\nr1 = 0.48r\norientation = out\n</plot>\n<plot>\nz = 2\nmax_gap = 0u\n')
            o.write('file = ' + p_gc_s + '\n')
            o.write(
                'color = green\nfill_color = green\nr0 = 0.25r\nr1 = 0.32r\norientation = out\n</plot>\n<plot>\nz = 2\nmax_gap = 0u\n')
            o.write('file = ' + n_gc_s + '\n')
            o.write(
                'color = orange\nfill_color = orange\nr0 = 0.18r\nr1 = 0.25r\norientation = out\n</plot>\n</plots>\n<<include etc/housekeeping.conf>>\n')
    else:
        if file:
            print(file)
        with open(output, 'w') as o:
            if int(type) == 1:
                multiplier = "1e-6"
                spacing1 = "20000b"
                spacing2 = "100000b"
            else:
                multiplier = "1e-3"
                spacing1 = "200b"
                spacing2 = "1000b"
                gsize = get_genome_size(kary)
                if gsize:
                    big, small = cal_space(gsize)
                    spacing2 = str(big) + 'b'
                    spacing1 = str(small) + 'b'

            o.write('<<include etc/colors_fonts_patterns.conf>>\n<image>\nangle_offset* = 98\n<<include etc/image.conf>>\n</image>\n\n')
            o.write('karyotype = ' + kary + '\n')
            o.write(
                'chromosomes_units = 100000\nchromosomes_display_default = yes\n\n<ideogram>\n<spacing>\n'
                'default = 0.005r\nbreak = 10r\naxis_break = yes\naxis_break_style = 2\naxis_break_at_edge = yes\n<break_style 1>\n'
                'stroke_color = black\nfill_color = blue\nthickness = 0.25r\nstroke_thickness = 2\n</break_style>\n<break_style 2>\n'
                'stroke_color = black\nstroke_thickness = 3\nthickness = 1.5r\n</break_style>\n'
                '</spacing>\nradius = 0.80r\nthickness = 1p\nfill = yes\nfill_color = deepskyblue\nstroke_color = black\n'  # thickness = 6p 改成 1p
                'stroke_thickness = 1p\nshow_label = no\nlabel_font = light\nlabel_radius = 1r + 110p\nlabel_size = 30\nlabel_parallel = no\n'
                '</ideogram>\n\nshow_ticks = yes\nshow_tick_labels = yes\n<ticks>\nskip_first_label = no\nskip_last_label = no\n'
                'radius = dims(ideogram,radius_outer)\ncolor = black\nthickness = 2p\nsize = 30p\nmultiplier = ' + multiplier + '\nformat = %.2f'
                '\n<tick>\nspacing = ' + spacing1 + '\nsize = 10p\nshow_label = no\nthickness = 3p\n</tick>\n<tick>\nspacing = ' + spacing2 + '\nsize = 20p\nshow_label = yes\n'
                                                                                                                                                                                                                                                              'label_size = 12p\nlabel_offset = 10p\nformat = %.2f\n</tick>\n</ticks>\n\n<highlights>\nz = 0\n<highlight>\n')
            # label_size 25 改成12
            o.write('file = ' + cog + '\n')
            o.write('r0 = 0.89r\nr1 = 0.99r\n</highlight>\n<highlight>\n')
            # o.write('file = ' + temp + '\n')
            # o.write('r0 = 0.89r\nr1 = 0.8901r\nfill_color = black\n</highlight>\n<highlight>\n')  #r1 = 0.891r  r1 = 0.8901r
            # o.write('r0 = 0.89r\nr1 = 0.891r\nstroke_color = black\nstroke_thickness = 1p\n</highlight>\n<highlight>\n')
            o.write('file = ' + an_cog + '\n')
            o.write('r0 = 0.79r\nr1 = 0.89r\n</highlight>\n')
            if file:
                file_list = file.split(',')
                div = (0.78 - 0.60) / len(file_list)
                n = 0
                for i in file_list:
                    o.write('<highlight>\nfile = ' + i + '\nr0 = ' + str(0.78 - n * div) + 'r\n' + 'r1 = ' + str(
                        0.78 - (n + 1) * div) + 'r\n</highlight>\n')
                    n = n + 1
            o.write('</highlights>\n\n')
            o.write('<plots>\ntype = line\nthickness = 1p\n<plot>\nz = 2\nmax_gap = 0u\n')
            o.write('file = ' + p_gc_c + '\n')
            o.write(
                'color = red\nfill_color = red\nr0 = 0.48r\nr1 = 0.58r\norientation = out\n</plot>\n<plot>\nz = 2\nmax_gap = 0u\n')
            o.write('file = ' + n_gc_c + '\n')
            o.write(
                'color = blue\nfill_color = blue\nr0 = 0.38r\nr1 = 0.48r\norientation = out\n</plot>\n<plot>\nz = 2\nmax_gap = 0u\n')
            o.write('file = ' + p_gc_s + '\n')
            o.write(
                'color = green\nfill_color = green\nr0 = 0.25r\nr1 = 0.32r\norientation = out\n</plot>\n<plot>\nz = 2\nmax_gap = 0u\n')
            o.write('file = ' + n_gc_s + '\n')
            o.write(
                'color = orange\nfill_color = orange\nr0 = 0.18r\nr1 = 0.25r\norientation = out\n</plot>\n</plots>\n<<include etc/housekeeping.conf>>\n')


def _main():
    parser = argparse.ArgumentParser(description='add gene_info in your table ')
    parser.add_argument('-k', '--kary', help="karyotype.txt")
    parser.add_argument('-c', '--cog', help="sense_strand_cog.txt")
    parser.add_argument('-t', '--temp', help="temp.txt")
    parser.add_argument('-ac', '--an_cog', help="antisense_strand_cog.txt")
    parser.add_argument('-pgc', '--p_gc_c', help="positive_gc_count.txt")
    parser.add_argument('-ngc', '--n_gc_c', help="negative_gc_count.txt")
    parser.add_argument('-pgs', '--p_gc_s', help="positive_gc_skew.txt")
    parser.add_argument('-ngs', '--n_gc_s', help="negative_gc_skew.txt")
    parser.add_argument('-f', '--file', default=None, help="other file")
    parser.add_argument('-ci', '--circle', default=None, help="circle type")
    parser.add_argument('-o', '--output', default=None, help="out file")
    parser.add_argument('-type', '--type', default=None, help="type 1 or 2 ")
    args = parser.parse_args()
    write_txt(args.output, args.kary, args.cog, args.temp, args.an_cog, args.p_gc_c, args.n_gc_c, args.p_gc_s,
              args.n_gc_s, args.file, args.type, args.circle)


if __name__ == "__main__":
    _main()
