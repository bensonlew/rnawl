#usr/bin/python

#zouguanqing 20190410

import svgwrite
import argparse
from compiler.ast import flatten
from biocluster.config import Config


def color_fun(labs):

    rna_colors = {"tRNA":"#FF0000","5S_rRNA":"#008000","16S_rRNA":"#0000FF","23S_rRNA":"#00BFFF"}
    def_circle = {"Def_circle":"#FF00FF"}
    island = {"Island":"#CDCD00"}
    prephage = {"prePhage":"#7B68EE"}
    isf = {"IS":"#1C1C1C"} #  28,28,28
    integron = {"integron":"#FF7F00"}  # 255,127,0

    ret =  {'rna':rna_colors, 'def_circle':def_circle, 'island':island, 'prephage':prephage, 'integron':integron, 'isf':isf}

    if 'COG' in labs:
        db_name = 'bacgenome'
        db = Config().get_mongo_client(mtype=db_name, ref=True)[Config().get_mongo_dbname(db_name, ref=True)]
        cog_color = {}
        for i in db.COG_color.find({'version':2}):
            func = db.eggNOG4_function_info.find_one({'function':i['fid']})
            key = i['fid'] + ':' + func['function_des']
            cog_color[key] = i['color']
        ret['COG'] = cog_color
    return ret


def draw_legend(lab2, colors,draw_class,start_x=100, start_y=100, square_size=15, font_size=8, interval_size=5,circle_id=None):
    position = [start_x, start_y]
    diff = square_size + interval_size
    current_circle = 0
    for n, labs in enumerate(lab2,1):
        print("labs: {}\n".format(labs))
        print("position: {}\n".format(position))
        # position = [position[0], position[1]+ diff]
        # position = [position[0], position[1]+ diff]
        # position_text = [position[0], position[1]+font_size]
        # for lab in labs:
        #     if lab not in colors.keys():
        #         raise(lab + ': name ERROR')
        #     if circle_id[lab] > current_circle:
        #         if n == 1:
        #             draw_class.add(draw_class.text('Circle 1-2:', position_text))
        #         else:
        #             draw_class.add(draw_class.text('Circle %s:'% (n+1), position_text))
        #
        #         current_circle = circle_id[lab]
        #         break

        # position = [position[0], position[1]+ diff]

        for lab in labs:
            print("lab: {}\n".format(lab))
            if lab not in colors.keys():
                raise(lab + ': name ERROR')
            for k in sorted(colors[lab].keys()):
                position = [position[0], position[1]+ diff]
                print("position2: {}\n".format(position))
                position_text = [position[0] + diff, position[1]+font_size]
                draw_class.add(draw_class.rect(position,(square_size,square_size),fill=colors[lab][k]))
                if k == 'prePhage':
                    k = 'prophage'
                draw_class.add(draw_class.text(k, position_text))


def main_fun(labs, out='legend.svg',start_x=100,start_y=100):
    lab2 = []
    circle_id = {}
    id = 1
    for i in labs:
        spi = i.split(',')
        for ei in spi:
            circle_id[ei] = id
        id +=1
        lab2.append(spi)
    merge_lab = flatten(lab2)
    colors = color_fun(merge_lab)
    legend=svgwrite.Drawing(filename=out,size=('100%','100%'),debug=True)
    draw_legend(lab2,colors,legend, start_x, start_y,circle_id=circle_id)
    legend.save()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in_str', help='value must in COG;def_circle,rna,island,prephage.Separate circles with semicolons.',required=True)
    parser.add_argument('-o','--out', help='out name',required=True)
    args = parser.parse_args()
    labs = args.in_str.split(';')
    #labs = ['def_circle','rna','island','prephage']
    main_fun(labs,args.out,1020,20)














