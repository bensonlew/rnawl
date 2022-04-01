# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from PIL import Image, ImageFont, ImageDraw
from Bio.Seq import Seq
import re


class MirandaPic(object):
    def __init__(self, infile, fonts):
        self.infile = infile
        self.fonts = fonts
        self.result = list()
        self.overall = dict()
        self.out = 'result_table.xls'

    def miranda_process(self):
        target = dict()
        with open(self.infile, 'r') as mi:
            line = mi.readline()
            while line:
                if line.strip().startswith("Query:"):
                    line2 = mi.next()
                    line3 = mi.next()
                    pairing = line + line2 + line3
                    q_seed = line.split()[2][-8:-1]
                    r_seed = line3.split()[2][-8:-1]
                    target['q_seed'] = q_seed[1:]
                    target['r_seed'] = r_seed[1:]
                    pattern = re.compile('(.*)' + q_seed[1:] + '(.*?)' + r_seed[1:] + '(.*)', re.S)
                    contents = pattern.search(pairing)
                    new_line = contents.group(2).split('\n')
                    target['pairing'] = [contents.group(1), new_line[0], new_line[1], new_line[2], contents.group(3)]
                    ref_1 = line3.split()[2][-1]
                    seed_seq = Seq(q_seed)
                    complement_seq = str(seed_seq.back_transcribe().complement())
                    mismatch = [x for x in range(len(r_seed)) if r_seed[x] != complement_seq[x]]
                    site_type = self.get_site_type(mismatch, ref_1.upper())
                    target['site_type'] = site_type
                    line = mi.next()
                    continue
                elif line.startswith(">>"):
                    info = line.strip().split()
                    self.overall = {'query': info[0].split('>>')[1],
                                    'ref': info[1],
                                    't_score': info[2],
                                    't_energy': info[3],
                                    'max_score': info[4],
                                    'max_energy': info[5],
                                    'query_len': info[7],
                                    'ref_len': info[8],
                                    'targets': info[9:],
                                    'mres': len(info[9:])}
                    break
                elif line.startswith(">"):
                    info = line.strip().split()
                    target['score'] = info[2]
                    target['energy'] = info[3]
                    target['start'] = info[6]
                    target['end'] = info[7]
                    self.result.append(target)
                    target = dict()
                    line = mi.next()
                    continue
                elif line.strip() == 'Scan Complete':
                    break
                else:
                    line = mi.next()
                    continue
        if len(self.result) > 0:
            self.output_generate()
            for each in self.result:
                self.draw_miranda_each(each)

    def get_site_type(self, mismatch, ref_1):
        num = len(mismatch)
        if ref_1 == 'A' and num == 0:
            site_type = '8mer site'
        elif ref_1 == 'A' and num == 1 and mismatch[0] == 0:
            site_type = '7mer-A1 site'
        elif ref_1 != 'A' and num == 0:
            site_type = '7mer-m8 site'
        elif ref_1 != 'A' and num == 1 and mismatch[0] == 0:
            site_type = '6mer site'
        elif ref_1 != 'A' and num == 1 and mismatch[0] == 6:
            site_type = 'Offset 6mer site'
        else:
            site_type = 'Unclassified'
        return site_type

    def output_generate(self):
        with open(self.out, 'w') as o:
            o.write('mirna\ttarget\ttotal_score\ttotal_energy\tmax_score\tmax_energy\tmi_length\ttarget_Length\tposition\n')
            o.write(self.overall['query'] + '\t' + self.overall['ref'] + '\t' + self.overall['t_score'] + '\t' +
                    self.overall['t_energy'] +'\t' + self.overall['max_score'] + '\t' + self.overall['max_energy'] +
                    '\t' + self.overall['query_len'] + '\t' + self.overall['ref_len'] + '\t' +
                    ';'.join(self.overall['targets']) + '\n')

    def draw_miranda_each(self, target):
        width = 2480
        height = 1500
        gap = 20
        height_rec = 80
        top = 30
        left = 60
        fonts = ImageFont.truetype(self.fonts, 60)
        out_pic = self.overall['query'] + '_vs_' + self.overall['ref'] + '_' + target['start'] + '.png'
        pics = Image.new("RGB", (width, height), (255, 255, 255))
        dr = ImageDraw.Draw(pics)
        # draw overview
        line = self.overall['query'] + '   vs.   ' + self.overall['ref']
        font_width, font_height = dr.textsize(line, fonts)
        dr.text(((width - font_width) / 2, top), line, fill="#000000", font=fonts)
        yaxis = top + font_height + gap
        dr.rectangle((left, yaxis, width-left, yaxis + height_rec), 'royalblue')
        dr.text(((width - dr.textsize('Overview', fonts)[0]) / 2, yaxis + (height_rec - font_height)/2), 'Overview', fill="white", font=fonts)
        yaxis += height_rec + 2*gap
        dr.text((150, yaxis), '5\'', fill='black', font=fonts)
        dr.text((2330, yaxis), '3\'', fill='black', font=fonts)
        yaxis += 2*font_height + gap
        font_width = dr.textsize('5\'', fonts)[0]
        dr.rectangle((150 + font_width, yaxis, 2330-font_width, yaxis+3), 'black')
        direction = 'up'
        for i in self.overall['targets']:
            top = 150 + font_width + int((2330-150-2*font_width) / float(self.overall['ref_len']) * float(i))
            label_x = top-dr.textsize(i, fonts)[0]/2
            if i == target['start']:
                colour = 'royalblue'
            else:
                colour = 'grey'
            if direction == 'up':
                dr.polygon((top, yaxis+3, top-24, yaxis+1728**0.5+3, top+24, yaxis+1728**0.5+3), fill=colour)  # equilateral triangle with side of 48
                dr.text((label_x, yaxis+1728**0.5+gap), i, fill=colour, font=fonts)
                direction = 'down'
            else:
                dr.polygon((top, yaxis, top - 24, yaxis - 1728 ** 0.5, top + 24, yaxis-1728 ** 0.5), fill=colour)
                dr.text((label_x, yaxis-1728**0.5-gap-font_height), i, fill=colour, font=fonts)
                direction = 'up'
        yaxis = int(yaxis + 1728**0.5 + 3 * gap + font_height)
        dr.text((left, yaxis), 'Query: ' + self.overall['query']+', Ref: '+self.overall['ref'], fill="#000000", font=fonts)
        dr.text((left, yaxis+font_height+gap), 'Total Score: '+self.overall['t_score']+', Total Energy: '+self.overall['t_energy'], fill="#000000", font=fonts)
        dr.text((left, yaxis + 2*font_height + 2*gap), 'Ref Length: ' + self.overall['ref_len'] + ', MREs: ' + str(self.overall['mres']), fill="#000000", font=fonts)
        yaxis += 3*font_height + 4*gap

        # draw target pair
        dr.rectangle((left, yaxis, width-left, yaxis+height_rec), 'royalblue')
        dr.text(((width - dr.textsize('2D Structure', fonts)[0]) / 2, yaxis+(height_rec-font_height)/2), '2D Structure', fill="white", font=fonts)
        yaxis += height_rec + 2*gap
        # fonts_seq = ImageFont.truetype("/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_042021/fonts/Courier-BOLD-3.ttf", 60)
        fonts_seq = fonts   # using the same font
        left1 = dr.textsize('   Query:    ', fonts_seq)[0]
        right1 = dr.textsize(target['pairing'][0] + target['q_seed'] + target['pairing'][1], fonts_seq)[0]
        dr.text((left+left1, yaxis), target['start'], fill='black', font=fonts)
        dr.text((left+right1, yaxis), target['end'], fill='black', font=fonts)
        xaxis = left + dr.textsize(target['pairing'][0], fonts_seq)[0]
        dr.text((xaxis+(dr.textsize(target['q_seed'], fonts_seq)[0] - dr.textsize('Seed', fonts)[0])/2, yaxis+0.5*gap), 'Seed', fill='red', font=fonts)
        yaxis += font_height + gap
        dr.text((left, yaxis), target['pairing'][0], fill='black', font=fonts_seq)
        dr.text((xaxis, yaxis), target['q_seed'], fill='red', font=fonts_seq)
        dr.text((xaxis, yaxis-0.5*font_height-2), '------', fill='red', font=fonts_seq)
        dr.text((xaxis + dr.textsize(target['q_seed'], fonts_seq)[0], yaxis), target['pairing'][1], fill='black', font=fonts_seq)
        dr.text((left, yaxis + font_height), target['pairing'][2], fill='black', font=fonts_seq)
        dr.text((left, yaxis + 2*font_height), target['pairing'][3], fill='black', font=fonts_seq)
        xaxis = left + dr.textsize(target['pairing'][3], fonts_seq)[0]
        dr.text((xaxis, yaxis + 2 * font_height), target['r_seed'], fill='royalblue', font=fonts_seq)
        dr.text((xaxis, yaxis + 2.5 * font_height-2), '------', fill='royalblue', font=fonts_seq)
        dr.text((xaxis+(dr.textsize(target['r_seed'], fonts_seq)[0] - dr.textsize(target['site_type'], fonts)[0])/2, yaxis + 3*font_height+ 0.5*gap), target['site_type'], fill='royalblue', font=fonts)
        dr.text((xaxis + dr.textsize(target['r_seed'], fonts_seq)[0], yaxis + 2 * font_height), target['pairing'][4], fill='black', font=fonts_seq)
        dr.text((left,  yaxis + 4*font_height + 2*gap), 'Score: ' + target['score']+', Energy: ' + target['energy'], fill='black', font=fonts)
        pics.save(out_pic)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='A script to visualise miRanda result.')
    parser.add_argument('-i', type=str, help='miRanda result file')
    parser.add_argument('-f', type=str, help='fonts file')
    args = parser.parse_args()
    mirandapic = MirandaPic(args.i, args.f)
    mirandapic.miranda_process()








