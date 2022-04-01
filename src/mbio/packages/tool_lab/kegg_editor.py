# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import os
import re
import gridfs
import subprocess
import sys
from biocluster.config import Config
from mbio.packages.ref_rna_v2.kegg_html import KeggHtml
from biocluster.config import Config
import json
import pandas as pd
from mbio.packages.rna.annot_config import AnnotConfig
from mbio.packages.ref_rna_v2.kegg_html import KeggHtml

class KeggEditor(object):
    def __init__(self, kegg_version=None):
        self.kegg_version = kegg_version
        # self.mongodb = Config().get_mongo_client(mtype='ref_rna', ref=True)[Config().get_mongo_dbname('ref_rna', ref=True)]
        # self.png_coll = self.mongodb.kegg_pathway_png_v1
        self.parafly = os.path.join(Config().SOFTWARE_DIR, 'bioinfo/denovo_rna_v2/trinityrnaseq-Trinity-v2.5.0/trinity-plugins/ParaFly-0.1.0/bin/ParaFly')
        self.png_modify_cmd = list()
        self.png2pdf = list()
        self.kegg_files_dict = AnnotConfig().get_file_dict(db="kegg", version=kegg_version)
        self.kegg_json = self.kegg_files_dict["br08901.json"]
        self.html_path = self.kegg_files_dict["html"]
        self.ko2gene = dict()

    def get_keggdb_paths(self):
        with open(self.kegg_json, 'rb') as f:
            root = json.load(f)
        classI = root['children']
        classII = []
        for i in classI:
            classII.extend(i['children'])
        classIII = []
        for i in classII:
            classIII.extend(i['children'])
        db_paths = ['map' + str(i['name']).split(' ')[0] for i in classIII]
        return db_paths

    def get_pic(self,r_path, map_path, path, kos_path, png_path, raw_html_path, png_bgcolor):
        map_id = re.sub("ko", "map", path)
        map_html = KeggHtml(version=self.kegg_version)
        map_html.color_bg[0] = png_bgcolor
        html_path = self.html_path + '/' + map_id + '.html'

        map_html.run(html_path, "pathways" + '/' + map_id + '.html', path + '.png', self.ko2gene)

        with open(raw_html_path + '/' + path + '.kgml', 'rb') as kgml:
            if len(kgml.readlines()) == 0:
                cmd = 'cp {}/{}.png {}'.format(raw_html_path, path, png_path)
            else:
                cmd = '{} {} {} {} {} {} {}'.format(
                    r_path,
                    map_path,
                    path,
                    kos_path,
                    png_path,
                    raw_html_path + '/' + path + '.kgml',
                    raw_html_path + '/' + path + '.png'
                )
        self.png_modify_cmd.append(cmd)

    def get_ko2gene(self, ko_gene_list):
        for k,g in ko_gene_list:
            if k in self.ko2gene:
                self.ko2gene[k] += "," + g
            else:
                self.ko2gene[k] = g
    def edit_kegg_pathway(self, all_html_path, r_path, map_path, image_magick, table1, table2, bg_color1, bg_color2, fg_colors1, fg_colors2, all_pathways="pathways"):
        if not os.path.exists(all_pathways):
            os.makedirs(all_pathways)

        f1 = {
            "up": fg_colors1[0],
            "down": fg_colors1[1]
        }

        if table1:
            df1 = pd.read_table(table1, header=0)
            df1["bg"] = bg_color1
            if "Regulate" in df1.columns:
                df1["fg"] = df1['Regulate'].map(lambda x:f1.get(x, ""))
            else:
                df1["fg"] = f1["up"]

            if "gene id" in df1.columns:
                self.get_ko2gene(zip(df1["KO id"], df1["gene id"]))

        # print table2
        if table2:
            print table2
            f2 = {
                "up": fg_colors2[0],
                "down": fg_colors2[1]
            }
            df2 = pd.read_table(table2, header=0)
            df2["bg"] = bg_color2
            if "Regulate" in df2.columns:
                df2["fg"] = df2['Regulate'].map(lambda x:f2.get(x, ""))
            else:
                df2["fg"] = f2["up"]
            print df2

            if "gene id" in df1.columns:
                self.get_ko2gene(zip(df2["KO id"], df2["gene id"]))

            df1 = df1.append(df2)

        # print df1
        # 去冗余合并颜色
        path_dict = dict()
        for line_dict in df1.to_dict("record"):
            if line_dict['Pathway id'] in path_dict:
                if line_dict["KO id"] in path_dict[line_dict['Pathway id']]:
                    ko_dict = path_dict[line_dict['Pathway id']][line_dict["KO id"]]
                    if line_dict['bg'] in ko_dict['bg']:
                        pass
                    else:
                        ko_dict['bg'].append(line_dict['bg'])

                    if line_dict['fg'] in ko_dict['fg']:
                        pass
                    else:
                        ko_dict['fg'].append(line_dict['fg'])
                else:
                    path_dict[line_dict['Pathway id']][line_dict["KO id"]] = {
                        "bg":[line_dict['bg']],
                        "fg":[line_dict['fg']],
                    }
            else:
                path_dict[line_dict['Pathway id']] = {
                    line_dict["KO id"]: {
                        "bg":[line_dict['bg']],
                        "fg":[line_dict['fg']],
                    }
                }

        for path in path_dict:
            with open(path + ".ko_color", 'w') as w:
                w.write('#KO\tbg\tfg\n')
                for k,v in path_dict[path].items():
                    w.write("\t".join([
                        k,
                        ",".join(v["bg"][:2]),
                        ",".join(v["fg"][:2]),
                    ]))
            png_path = all_pathways + '/' + path + '.png'
            pdf_path = all_pathways + '/' + path + '.pdf'
            html_path = all_pathways + '/' + path + '.html'
            png_bgcolor = ",".join(v["bg"][:2])
            self.get_pic(r_path, map_path, path, path + ".ko_color" , png_path, all_html_path, png_bgcolor)

        with open('cmd_getpic.sh', 'wb') as f:
            f.write('\n'.join(self.png_modify_cmd) + "\n")
        cmd = '{} -c {} -CPU {} -v -shuffle'.format(self.parafly, 'cmd_getpic.sh', 10)
        print cmd
        os.system(cmd)
        # with open(prefix + 'cmd_png2pdf.sh', 'wb') as f:
        #     f.write('\n'.join(self.png2pdf))
        # cmd = '{} -c {} -CPU {} -v -shuffle'.format(self.parafly, prefix +  'cmd_png2pdf.sh', 10)
        # print cmd
        # os.system(cmd)

if __name__ == '__main__':
    if sys.argv[6] == "None":
        table2 = None
    else:
        table2 = sys.argv[6]
    KeggEditor(kegg_version="202007").edit_kegg_pathway(
        all_html_path=sys.argv[1],
        r_path=sys.argv[2],
        map_path=sys.argv[3],
        image_magick=sys.argv[4],
        table1=sys.argv[5],
        table2=table2,
        bg_color1=sys.argv[7],
        bg_color2=sys.argv[8],
        fg_colors1=sys.argv[9].split(","),
        fg_colors2=sys.argv[10].split(",")
    )
