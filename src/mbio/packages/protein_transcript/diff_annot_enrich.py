#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author fengyitong 2019-02

import sys
sys.path.insert(0, '/mnt/ilustre/users/yitong.feng/scripts/run_commands')
from controller import Controller
sys.path.insert(0, '/mnt/ilustre/users/hui.wan/.local/lib/python2.7/site-packages/Mako-1.0.6-py2.7.egg')
from mako.template import Template
import os
import argparse
import glob


class DiffAnnotEnrich(Controller):
    def __init__(self, go_list, go_obo, go_level, path_txt, kegg_dir, orgc, pathway_table, diff_path, out):
        super(DiffAnnotEnrich, self).__init__()
        self.go_list = os.path.abspath(go_list)
        if not os.path.exists(self.go_list):
            super(DiffAnnotEnrich, self).end(normal=1, out='你传入的go list文件不存在')
        with open(self.go_list, 'r') as gor:
            self.acc2go = {line.split('\t')[0]:line.strip() for line in gor if line.strip()}
        # 确定画炫图时候图片的大小，避免accession号过长而显示不全的问题
        self.acc_tolong = False
        for acc in self.acc2go:
            if len(acc) > 15:
                self.acc_tolong = True
        self.go_obo = os.path.abspath(go_obo)
        if not os.path.exists(self.go_obo):
            super(DiffAnnotEnrich, self).end(normal=1, out='你传入的go_obo文件不存在')
        self.go_level = os.path.abspath(go_level)
        if not os.path.exists(self.go_level):
            super(DiffAnnotEnrich, self).end(normal=1, out='你传入的go_level文件不存在')
        self.path_txt = os.path.abspath(path_txt)
        if not os.path.exists(self.path_txt):
            super(DiffAnnotEnrich, self).end(normal=1, out='你传入的path_txt文件不存在')
        with open(self.path_txt, 'r') as keggr:
            self.acc2kegg = {line.split('\t')[0]:line.strip() for line in keggr if line.strip()}
        self.kegg_dir = os.path.abspath(kegg_dir)
        if not os.path.exists(self.kegg_dir):
            super(DiffAnnotEnrich, self).end(normal=1, out='你传入的kegg_dir文件不存在')
        self.pathway_table = os.path.abspath(pathway_table)
        if not os.path.exists(self.pathway_table):
            super(DiffAnnotEnrich, self).end(normal=1, out='你传入的pathway_table文件不存在')
        self.out = os.path.abspath(out)
        if not os.path.exists(self.out):
            try:
                os.mkdir(self.out)
            except:
                super(DiffAnnotEnrich, self).end(normal=1, out='输出文件夹路径不对')
        path_xls = os.path.join(self.kegg_dir, orgc + '.pathway.xls')
        if os.path.exists(path_xls):
            self.typeii = os.path.join(self.out, 'pathwayId_pathwayName_typeII_typeI')
            cmd = """awk 'BEGIN{FS=OFS="\t"}{if (NR==1){print "pathwayId","#Term","typeII","typeI"}else{print $3,$4,$2,$1}}' ${path} > ${type_}"""
            # cmd = r"""awk "BEGIN{FS=OFS="\t"}{if (NR==1){print "pathwayId","#Term","typeII","typeI"}else{print $3,$4,$2,$1}}" ${path} > ${type_}"""
            cmd = Template(cmd)
            cmd = cmd.render(path=path_xls,
                             type_=self.typeii
                             )
            os.system(cmd)
        else:
            self.typeii = '/mnt/ilustre/users/yitong.feng/scripts/diff_annot_enrich/pathwayId_pathwayName_typeII_typeI'
        self.diff_path = os.path.abspath(diff_path)
        if not os.path.exists(self.diff_path):
            super(DiffAnnotEnrich, self).end(normal=1, out='你传入的diff结果路径不存在')
        self.DElists = glob.glob(self.diff_path + '/*.DE.list')
        if not self.DElists:
            super(DiffAnnotEnrich, self).end(normal=1, out='你传入的diff结果路径下没有DE文件')
        self.kegg_list = list()
        self.go_list = list()
        # self.cmp_kegg = list()
        self.cmp_go = list()
        # self.cmp2kegg = dict()
        self.cmp2go = dict()

    def extract(self, infile, outfile, type_ = 'go'):
        with open(infile, 'r') as r, open(outfile, 'w') as w:
            for line in r:
                if type_ == 'go':
                    try:
                        w.write(self.acc2go[line.strip()] + '\n')
                    except:
                        pass
                else:
                    try:
                        w.write(self.acc2kegg[line.strip()] + '\n')
                    except:
                        pass

    def add_commands(self):
        for delist in self.DElists:
            cmp = os.path.basename(delist).split('.DE.list')[0]
            self.extract(os.path.join(self.diff_path, cmp + '.DE.list'), os.path.join(self.out, cmp + '.DE.GO.list'), type_='go')
            self.extract(os.path.join(self.diff_path, cmp + '.up.list'), os.path.join(self.out, cmp + '.up.GO.list'), type_='go')
            self.extract(os.path.join(self.diff_path, cmp + '.down.list'), os.path.join(self.out, cmp + '.down.GO.list'), type_='go')
            self.extract(os.path.join(self.diff_path, cmp + '.DE.list'), os.path.join(self.out, cmp + '.DE.pathway.txt'), type_='kegg')
            self.extract(os.path.join(self.diff_path, cmp + '.up.list'), os.path.join(self.out, cmp + '.up.pathway.txt'), type_='kegg')
            self.extract(os.path.join(self.diff_path, cmp + '.down.list'), os.path.join(self.out, cmp + '.down.pathway.txt'), type_='kegg')
            go_de = self.add_command(name='go_+_DE_+_' + cmp)
            self.go_list.append(go_de)
            go_up = self.add_command(name='go_+_up_+_' + cmp)
            self.go_list.append(go_up)
            go_down = self.add_command(name='go_+_down_+_' + cmp)
            self.go_list.append(go_down)
            kegg_de = self.add_command(name='kegg_+_DE_+_' + cmp)
            self.kegg_list.append(kegg_de)
            kegg_up = self.add_command(name='kegg_+_up_+_' + cmp)
            self.kegg_list.append(kegg_up)
            kegg_down = self.add_command(name='kegg_+_down_+_' + cmp)
            self.kegg_list.append(kegg_down)
            # cmp_kegg = self.add_command(name='kegg_' + cmp)
            # self.cmp_kegg.append(cmp_kegg)
            cmp_go = self.add_command(name='go_' + cmp)
            self.cmp_go.append(cmp_go)
            # self.cmp2kegg[cmp] = dict(
            #     up=kegg_up,
            #     down= kegg_down,
            #     de=kegg_de
            # )
            self.cmp2go[cmp] = dict(
                up=go_up,
                down=go_down,
                de=go_de,
                cmp_go=cmp_go
            )

    def how_run(self):
        self.after_some(self.go_list, self.run_cmp_go)
        end_list = self.cmp_go + self.kegg_list
        self.after_some(end_list, self.set_output)

    def end(self):
        super(DiffAnnotEnrich, self).end(out='your program finished')

    def run_kegg(self):
        wp = 25
        hp = 27
        if self.acc_tolong:
            wp = 60
            hp = 62
        for kegg in self.kegg_list:
            _, regu, cmp = kegg.name.split('_+_')
            try:
                # os.link(os.path.join(self.diff_path, '%s.diff.exp.xls'%cmp), os.path.join(kegg.work_dir, '%s.diff.exp.xls'%cmp))
                cmd = 'cp ' + os.path.join(self.diff_path, '%s.diff.exp.xls'%cmp) + ' ' + os.path.join(kegg.work_dir, '%s.diff.exp.xls'%cmp)
                os.system(cmd)
            except:
                pass
            if regu.lower() == 'up':
                cmd = r"""#${bin}/KEGG_all_annotation_col.py -paths ${KEGG_dir} -pathwaytxt ${pathwaytxt} -out ${cmp}.${regu}.paths
touch temp.list
${bin}/KEGG_diff_annotation_col.py -paths ${KEGG_dir} -pathwaytxt ${pathwaytxt} -diffxls ${diff_path}/${cmp}.diff.exp.xls -uplist ${diff_path}/${cmp}.up.list -downlist temp.list -out ${cmp}.up.paths
${bin}/KEGG_enrichment_xls.py -all ${pathway_table} -diff ${cmp}.${regu}.paths/pathway_table.xls -p_class ${typeiii} -out_xls ${cmp}.${regu}.kegg_enrichment.xls
awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$6,$8}' ${cmp}.${regu}.kegg_enrichment.xls > plot.${cmp}.${regu}.kegg_enrichment.xls
${bin}/plot_enrichment_bar.py -i plot.${cmp}.${regu}.kegg_enrichment.xls -o ${cmp}.${regu}.kegg_enrichment.pdf
${bin}/kegg_enrichment_bubble.py ${cmp}.${regu}.kegg_enrichment.xls ${cmp}.${regu}.kegg_enrichment_bubble.pdf
${bin}/keggplot_circle.py -diffxls ${cmp}.diff.exp.xls -diff_enrich ${cmp}.${regu}.kegg_enrichment.xls -out ${cmp}.${regu}.KEGGplot.pdf -wp ${wp} -hp ${hp}
"""
            elif regu.lower() == 'down':
                cmd = r"""#${bin}/KEGG_all_annotation_col.py -paths ${KEGG_dir} -pathwaytxt ${pathwaytxt} -out ${cmp}.${regu}.paths
touch temp.list
${bin}/KEGG_diff_annotation_col.py -paths ${KEGG_dir} -pathwaytxt ${pathwaytxt} -diffxls ${diff_path}/${cmp}.diff.exp.xls -uplist temp.list -downlist ${diff_path}/${cmp}.down.list -out ${cmp}.down.paths
${bin}/KEGG_enrichment_xls.py -all ${pathway_table} -diff ${cmp}.${regu}.paths/pathway_table.xls -p_class ${typeiii} -out_xls ${cmp}.${regu}.kegg_enrichment.xls
awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$6,$8}' ${cmp}.${regu}.kegg_enrichment.xls > plot.${cmp}.${regu}.kegg_enrichment.xls
${bin}/plot_enrichment_bar.py -i plot.${cmp}.${regu}.kegg_enrichment.xls -o ${cmp}.${regu}.kegg_enrichment.pdf
${bin}/kegg_enrichment_bubble.py ${cmp}.${regu}.kegg_enrichment.xls ${cmp}.${regu}.kegg_enrichment_bubble.pdf
${bin}/keggplot_circle.py -diffxls ${cmp}.diff.exp.xls -diff_enrich ${cmp}.${regu}.kegg_enrichment.xls -out ${cmp}.${regu}.KEGGplot.pdf -wp ${wp} -hp ${hp}
"""
            else:
                regu = 'DE'
                try:
                    os.link(self.path_txt, os.path.join(kegg.work_dir, 'pathway.txt'))
                except:
                    pass
                cmd = r"""${bin}/KEGG_diff_annotation_col.py -paths ${KEGG_dir} -pathwaytxt pathway.txt -diffxls ${diff_path}/${cmp}.diff.exp.xls -uplist ${diff_path}/${cmp}.up.list -downlist ${diff_path}/${cmp}.down.list -out ${cmp}.DE.paths
${bin}/KEGG_enrichment_xls.py -all ${pathway_table} -diff ${cmp}.${regu}.paths/pathway_table.xls -p_class ${typeiii} -out_xls ${cmp}.${regu}.kegg_enrichment.xls
awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$6,$8}' ${cmp}.${regu}.kegg_enrichment.xls > plot.${cmp}.${regu}.kegg_enrichment.xls
${bin}/plot_enrichment_bar.py -i plot.${cmp}.${regu}.kegg_enrichment.xls -o ${cmp}.${regu}.kegg_enrichment.pdf
${bin}/kegg_enrichment_bubble.py ${cmp}.${regu}.kegg_enrichment.xls ${cmp}.${regu}.kegg_enrichment_bubble.pdf
${bin}/keggplot_circle.py -diffxls ${cmp}.diff.exp.xls -diff_enrich ${cmp}.${regu}.kegg_enrichment.xls -out ${cmp}.${regu}.KEGGplot.pdf -wp ${wp} -hp ${hp}
"""
            cmd = Template(cmd)
            cmd = cmd.render(bin='/mnt/ilustre/users/ting.kuang/ITRAQ/bin',
                            KEGG_dir=self.kegg_dir,
                            pathwaytxt=os.path.join(self.out, cmp + '.%s.pathway.txt'%regu),
                             cmp=cmp,
                             regu=regu,
                             typeiii=self.typeii,
                             diff_path=self.diff_path,
                             pathway_table=self.pathway_table,
                             wp = wp,
                             hp = hp
                           )
            params = dict(
                cmd=cmd,
                node=4,
                memory=20
            )
            kegg.set_params(params)
            kegg.run()

    def run_go(self):
        wp = 25
        hp = 27
        if self.acc_tolong:
            wp = 60
            hp = 62
        for go in self.go_list:
            _, regu, cmp = go.name.split('_+_')
            if regu == 'de':
                regu = 'DE'
            try:
                # os.link(os.path.join(self.diff_path, '%s.diff.exp.xls'%cmp), os.path.join(go.work_dir, '%s.diff.exp.xls'%cmp))
                cmd = 'cp ' + os.path.join(self.diff_path, '%s.diff.exp.xls' % cmp) + ' ' + os.path.join(go.work_dir,
                                                                                                         '%s.diff.exp.xls' % cmp)
                os.system(cmd)
            except:
                pass
            cmd = r"""${bin}/go_annotation.py -go_obo ${go_obo} -i ${out}/${cmp}.${regu}.GO.list -o ${cmp}.${regu}.GO.level.xls
${bin}/GO_enrichment_xls.py -all ${go_level} -diff ${cmp}.${regu}.GO.level.xls -out_xls ${cmp}.${regu}.go_enrichment.xls
awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$6,$8}' ${cmp}.${regu}.go_enrichment.xls > plot.${cmp}.${regu}.go_enrichment.xls
${bin}/plot_enrichment_bar.py -i plot.${cmp}.${regu}.go_enrichment.xls -o ${cmp}.${regu}.go_enrichment.pdf
/mnt/ilustre/users/yitong.feng/scripts/diff_annot_enrich/go_enrichment_bubble.py ${cmp}.${regu}.go_enrichment.xls ${cmp}.${regu}.go_enrichment_bubble.pdf
/mnt/ilustre/users/ting.kuang/ALL-SCRIPT/goplot_circle.py -diffxls ${cmp}.diff.exp.xls -diff_enrich ${cmp}.${regu}.go_enrichment.xls -out ${cmp}.${regu}.GOplot.pdf -wp ${wp} -hp ${hp}
"""
            cmd = Template(cmd)
            cmd = cmd.render(bin='/mnt/ilustre/users/ting.kuang/ITRAQ/bin',
                             go_obo=self.go_obo,
                             out=self.out,
                             cmp=cmp,
                             regu=regu,
                             go_level=self.go_level,
                             wp = wp,
                             hp = hp
                             )
            params = dict(
                cmd=cmd,
                node=4,
                memory=20
            )
            go.set_params(params)
            go.run()

    def run_cmp_go(self):
        for cmp in self.cmp2go:
            try:
                os.link(os.path.join(self.cmp2go[cmp]['up'].work_dir, '%s.up.GO.level.xls'%(cmp)), os.path.join(self.cmp2go[cmp]['cmp_go'].work_dir, '%s.up.GO.level.xls'%(cmp)))
                os.link(os.path.join(self.cmp2go[cmp]['down'].work_dir, '%s.down.GO.level.xls'%(cmp)), os.path.join(self.cmp2go[cmp]['cmp_go'].work_dir, '%s.down.GO.level.xls'%(cmp)))
            except:
                pass
            cmd = '/mnt/ilustre/users/ting.kuang/ITRAQ/bin/go_level_up_down.pl -u {cmp}.up.GO.level.xls -d {cmp}.down.GO.level.xls -lg 2 -o {cmp}.GO.level2.up-down.pdf'.format(cmp=cmp)
            params = dict(
                cmd=cmd,
                node=1,
                memory=1
            )
            self.cmp2go[cmp]['cmp_go'].set_params(params)
            self.cmp2go[cmp]['cmp_go'].run()

    def set_output(self):
        all_c = self.go_list + self.kegg_list + self.cmp_go
        for c in all_c:
            for file in os.listdir(c.work_dir):
                if (file.endswith('.pdf') or file.endswith('.xls')) and not file.startswith('tmp'):
                    try:
                        os.link(os.path.join(c.work_dir, file), os.path.join(self.out, file))
                    except:
                        pass
                if os.path.isdir(os.path.join(c.work_dir, file)):
                    try:
                        os.system('cp -r %s %s' %(os.path.join(c.work_dir, file), self.out))
                    except:
                        pass
        self.end()

    def run(self):
        self.add_commands()
        self.how_run()
        self.run_kegg()
        self.run_go()
        self.fire()

if __name__ == '__main__':
    # def __init__(self, go_list, go_obo, go_level, path_txt, kegg_dir, orgc, pathway_table, diff_path, out):
    parser = argparse.ArgumentParser(description="The script to run diff annot and enrich for itraq")
    parser.add_argument("-go_list", type=str, required=True, help="the go_list file")
    parser.add_argument("-go_obo", type=str, required=True, help="the go_obo database")
    parser.add_argument("-go_level", type=str, required=True, help="the go annot result")
    parser.add_argument("-path_txt", type=str, required=True, help="the path_txt file")
    parser.add_argument("-kegg_dir", type=str, required=True, help="the kegg pictures")
    parser.add_argument("-pathway_table", type=str, required=True, help="the kegg annot result")
    parser.add_argument("-orgc", type=str, default='kegg', help='the organsim of the species')
    parser.add_argument("-diff_path", type=str, required=True, help="the diff result path")
    parser.add_argument("-out", type=str, default=os.path.join(os.getcwd(), 'diff_annot_enrich_results'))

    args = parser.parse_args()
    diffae = DiffAnnotEnrich(args.go_list, args.go_obo, args.go_level, args.path_txt, args.kegg_dir, args.orgc, args.pathway_table, args.diff_path, args.out)
    diffae.run()