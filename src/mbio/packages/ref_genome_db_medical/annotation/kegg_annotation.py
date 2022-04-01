# -*- coding: utf-8 -*-
# __author__ = 'zengjing'

from biocluster.config import Config
import xml.etree.ElementTree as ET
import re
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis_t import KGMLCanvas
from mbio.packages.ref_rna_v2.kegg_html import KeggHtml
from reportlab.lib import colors
import collections
import json
from itertools import islice
import subprocess
import gridfs
import os
import sys

class KeggAnnotation(object):
    def __init__(self):
        '''
        设置数据库，连接到mongod数据库，涉及kegg_ko，kegg_gene，kegg_pathway_png三个collections
        '''
        self.client = Config().get_mongo_client(mtype='ref_rna', ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname('ref_rna', ref=True)]
        self.gene_coll = self.mongodb.kegg_gene_v1
        self.ko_coll = self.mongodb.kegg_ko_v1
        self.png_coll = self.mongodb.kegg_pathway_png_v1
        self.path = collections.defaultdict(str)
        self.ko2gene = dict()
        self.kegg_json = Config().SOFTWARE_DIR + '/database/KEGG/br08901.json'
        self.kegg_pathway = Config().SOFTWARE_DIR + '/database/KEGG/pathway'
        self.kegg_ko = Config().SOFTWARE_DIR + '/database/KEGG/ko_des'
        self.gene_ko = Config().SOFTWARE_DIR + '/database/KEGG/ko_genes.list'
        self.ko2name = dict()
        self.gloabl = [
            'map01100', 'map01110', 'map01120', 'map01130', 'map01200', 'map01210', 'map01212', 'map01230', 'map01220'
        ]

    def get_ko2path(self):
        '''
        获取ko与pathway对应关系
        '''
        ko2path = dict()
        for line in open(self.kegg_pathway):
            cols = line.strip().split('\t')
            path = cols[0].split(':')[1]
            ko = cols[1].split(':')[1]
            if path in self.gloabl:
                continue
            if ko in ko2path:
                if path in ko2path[ko]:
                    pass
                else:
                    ko2path[ko].add(path)
            else:
                ko2path[ko] = set([path])
        return ko2path

    def get_ko2name(self):
        '''
        获取ko与描述对应关系
        '''
        ko2des = dict()
        for line in open(self.kegg_ko):
            cols = line.strip().split('\t')
            ko = cols[0].split(':')[1]
            name = cols[1].split(';')[0]
            if ko in ko2des:
                pass
            else:
                ko2des[ko] = name
        self.ko2name = ko2des
        return ko2des

    def get_gene2ko(self):
        '''
        获取gene与ko对应关系
        '''
        gene2ko = dict()
        for line in open(self.gene_ko, 'r'):
            cols = line.strip().split("\t")
            ko = cols[0].split(":")[1]
            gene = cols[1]
            if gene in gene2ko:
                gene2ko[gene].add(ko)
            else:
                gene2ko[gene] = set([ko])
        return gene2ko

    def get_kegg_class(self):
        '''
        获取kegg分类注释
        '''
        map2class = dict()
        with open(self.kegg_json) as f:
            root = json.load(f)
        classI = root['children']
        for classI_child in classI:
            classI_name = classI_child['name']
            for classII_child in classI_child['children']:
                classII_name = classII_child['name']
                for classIII_child in classII_child['children']:
                    classIII_name = classIII_child['name']
                    names = classIII_name.split('  ')
                    path_id = 'map' + names[0]
                    name = names[1]
                    map2class[path_id] = [classI_name, classII_name, name]
        return map2class

    def pathSearch(self, blast_xml, kegg_table, taxonomy=None):
        '''
        输入blast.xml
        输出kegg_table.xls
        '''
        ko2path = self.get_ko2path()
        gene2ko = self.get_gene2ko()
        ko2name = self.get_ko2name()
        path_db = self.get_keggdb_paths()
        map2class = self.get_kegg_class()
        ko_list = list()
        if taxonomy:
            for line in open(taxonomy):
                line = line.strip()
                ko_list.append(line)
        tablefile = open(kegg_table, 'wb')
        tablefile.write('#Query\tKO_ID (Gene id)\tKO_name (Gene name)\tHyperlink\tPaths\n')
        docment = ET.parse(blast_xml)
        root = docment.getroot()

        iterns = root.find('BlastOutput_iterations')
        for itern in iterns:
            query = itern.find('Iteration_query-def').text.split()[0]
            iter_hits = itern.find('Iteration_hits')
            hits = iter_hits.findall('Hit')
            if len(hits) > 0:
                mark = 0
                ko, ko_name, ko_hlink, path_list = list(), list(), list(), list()
                for hit in hits:
                    mark += 1
                    if mark == 6:
                        break
                    gid = hit.find('Hit_id').text
                    if not gid in gene2ko:
                        continue
                    koids = gene2ko[gid]
                    for item in koids:
                        ko.append(item)
                        if item in ko2path:
                            pids = ko2path[item]
                            for index, i in enumerate(pids):
                                map_id = re.sub('ko', 'map', i)
                                if map_id not in self.gloabl:
                                    self.path[map_id] = map2class[map_id]
                                    path_list.append(map_id)
                ko_list = list(set(ko))
                ko = ';'.join(ko_list)
                for i in ko.split(';'):
                    if self.ko2gene.has_key(i):
                        self.ko2gene[i] += '|' + query
                    else:
                        self.ko2gene[i] = 'accession: ' + query
                ko_name = ';'.join([ko2name[x] if x in ko2name else '' for x in ko_list])
                ko_hlink = 'http://www.genome.jp/dbget-bin/www_bget?ko:' + ko
                path = ';'.join(list(set(path_list)))
                if not path:
                    path = '\t'
                if ko:
                    if ko_name:
                        tablefile.write(query + '\t' + ko + '\t' + ko_name + '\t' + ko_hlink + '\t' + path + '\n')
                    else:
                        print '没有在kegg_ko数据库找到{}'.format(ko_name)
                else:
                    print '{}没有在数据库kegg_gene找到相应的koid'.format(gid)
            else:
                print '没有找到在该query下对应的基因信息'
        root.clear()
        print 'pathSearch finished'

    def merge_known(self, known_ko_annot, kegg_table):
        '''
        合并已知注释
        '''
        known_ko = dict()
        annot_ko = dict()
        self.ko2gene = dict()
        for line in open(known_ko_annot):
            cols = line.strip().split('\t')
            known_ko[cols[1]] = line.strip()

        header = str()
        with open(kegg_table) as f_annot:
            header = f_annot.read()
            for line in f_annot:
                cols = line.strip().split('\t')
                annot_ko[cols[0]] = line.strip()

        with open(kegg_table, 'w') as w_annot:
            w_annot.write(header)
            for tran_id in list(set(known_ko.keys() + annot_ko.keys())):
                if tran_id in known_ko.keys():
                    ko_annot = known_ko[tran_id].split('\t')
                    ko = ko_annot[3]
                    for i in ko.split(';'):
                        if self.ko2gene.has_key(i):
                            self.ko2gene[i] += '|' + tran_id
                        else:
                            self.ko2gene[i] = 'accession: ' + tran_id

                    ko_hlink = 'http://www.genome.jp/dbget-bin/www_bget?ko:' + ko
                    if ko in self.ko2name:
                        name = self.ko2name[ko]
                    else:
                        name = str()
                        print '{} do not have name'.format(ko)
                    if len(ko_annot) >= 5:
                        maps = ';'.join(['map' + x for x in ko_annot[4].split(';') if 'map' + x not in self.gloabl])
                    else:
                        maps = str()
                    w_annot.write('\t'.join([
                        tran_id,
                        ko,
                        name,
                        ko_hlink,
                        maps
                    ]) + '\n')
                else:
                    w_annot.write(annot_ko[tran_id] + '\n')
                    ko = annot_ko[tran_id].split('\t')[2]
                    for i in ko.split(';'):
                        if self.ko2gene.has_key(i):
                            self.ko2gene[i] += '|' + tran_id
                        else:
                            self.ko2gene[i] = 'accession: ' + tran_id

    def pathSearch_upload(self, kegg_ids, kegg_table, taxonomy=None):
        '''
        输入基因或转录本id对应的K编号文件kegg.list
        输出kegg_table.xls
        '''
        tablefile = open(kegg_table, 'wb')
        ko_list = list()
        if taxonomy:
            for line in open(taxonomy):
                line = line.strip()
                ko_list.append(line)
        tablefile.write('#Query\tKO_ID (Gene id)\tKO_name (Gene name)\tHyperlink\tPaths\n')
        for line in open(kegg_ids):
            ko, ko_name, ko_hlink, path = list(), list(), list(), list()
            line = line.strip().split('\t')
            query = line[0]
            kos = line[1].split(';')
            for ko_id in kos:
                ko.append(ko_id)
                result = self.ko_coll.find_one({'ko_id': ko_id})
                if result:
                    ko_name.append(result['ko_name'])
                    pids = result['pathway_id']
                    for index, i in enumerate(pids):
                        if ko_list:
                            if i in ko_list:
                                map_id = re.sub('ko', 'map', i)
                                if map_id not in self.gloabl:
                                    self.path[map_id] = result['pathway_category'][index]
                                    path.append(map_id)
                        else:
                            map_id = re.sub('ko', 'map', i)
                            if map_id not in self.gloabl:
                                self.path[map_id] = result['pathway_category'][index]
                                path.append(map_id)
                else:
                    print '没有在kegg_ko数据库找到{}'.format(ko_id)
            ko = ';'.join(ko)
            ko_name = ';'.join(ko_name)
            ko_hlink = 'http://www.genome.jp/dbget-bin/www_bget?ko:' + ko
            path = ';'.join(path)
            if not path:
                path = '\t'
            if ko:
                if ko_name:
                    tablefile.write(query + '\t' + ko + '\t' + ko_name + '\t' + ko_hlink + '\t' + path + '\n')
                else:
                    print '没有在kegg_ko数据库找到{}'.format(ko_id)
        print 'pathSearch finished'

    def get_keggdb_paths(self):
        with open(self.kegg_json) as f:
            root = json.load(f)
        classI = root['children']
        classII = list()
        for i in classI:
            classII.extend(i['children'])
        classIII = list()
        for i in classII:
            classIII.extend(i['children'])
        db_paths = ['map' + str(i['name']).split(' ')[0] for i in classIII]
        return db_paths

    def pathTable(self, r_path, map_path, kegg_table, pathway_path, pidpath, link_bgcolor, png_bgcolor, pathwaydir, image_magick):
        '''
        根据pathSearch生成的kegg_table.xls统计pathway的信息
        输入文件为kegg_table.xls
        输出文件为pathway_table.xls,pid.txt
        '''
        if not os.path.exists(pathwaydir):
            os.makedirs(pathwaydir)
        path_table_xls = open(pathway_path, 'wb')
        pid_txt = open(pidpath, 'wb')
        header_line = 'Pathway' + '\t' + 'First Category' + '\t' + 'Second Category' + '\t' + 'Pathway_definition' + '\t' + 'num_of_seqs' + '\t' + 'seqs_kos/gene_list' + '\t' + 'pathway_imagename' + '\t' + 'Hyperlink' + '\n'
        path_table_xls.write(header_line)
        path_table = collections.defaultdict(list)
        kegg_table = islice(open(kegg_table), 1, None)
        kegg = [i.strip('\n').split('\t') for i in kegg_table]
        table = [(i[0] + '(' + i[1] + ')', i[4]) for i in kegg]
        for i in table:
            if i[1] == str():
                continue
            for path in i[1].split(';'):
                print 'path is {}'.format(path)
                path_table[path].append(i[0])
        path_db = self.get_keggdb_paths()
        sorted_paths = sorted(path_table.keys(), key=lambda x: path_db.index(x))
        for key in sorted_paths:
            if key:
                pid = re.sub('map', 'ko', key)
                if key in self.path:
                    definition = self.path[key][2]
                else:
                    definition = str()
                koids = list()
                for i in path_table[key]:
                    for j in i.split('(')[1][0:-1].split(';'):
                        koids.append(j)
                koids = set(koids)
                koid_str = ';'.join(koids)
                ko_color = list()
                fgcolor = 'NA'
                kos_path = os.path.join(os.getcwd(), 'KOs.txt')
                with open(kos_path, 'w') as w:
                    w.write('#KO\tbg\tfg\n')
                    for k in koids:
                        ko_color.append(k + '%09' + link_bgcolor)
                        w.write(k + '\t' + png_bgcolor + '\t' + fgcolor + '\n')
                png_path = pathwaydir + '/' + key + '.png'
                pdf_path = pathwaydir + '/' + key + '.pdf'
                html_path = pathwaydir + '/' + key + '.html'
                self.get_pic(r_path, map_path, key, kos_path, png_path, html_path, png_bgcolor)
                if image_magick:
                    cmd = '{} -flatten -quality 100 -density 130 -background white {} {}'.format(
                        image_magick, png_path, pdf_path
                    )
                    try:
                        subprocess.check_output(cmd, shell=True)
                    except subprocess.CalledProcessError:
                        print '图片格式pdf转png出错'
                link = 'http://www.genome.jp/dbget-bin/show_pathway?' + key + '/' + '/'.join(ko_color)
                pid_txt.write(key + '\t' + koid_str + '\n')
                result = self.ko_coll.find_one({'pathway_id': {'$in': [pid]}})
                if result:
                    pids = result['pathway_id']
                    for index, i in enumerate(pids):
                        if i == pid:
                            category = result['pathway_category'][index]
                            layer_1st = category[0]
                            layer_2nd = category[1]
                    num_of_seqs = len(path_table[key])
                    geneids = [j.split('(')[0] for j in path_table[key]]
                    genes = ';'.join(geneids)
                    path_image = key + '.png'
                    line = key + '\t' + layer_1st + '\t' + layer_2nd + '\t' + definition + '\t' + str(num_of_seqs) + '\t' + genes + '\t' + path_image + '\t' + link + '\n'
                    path_table_xls.write(line)
            else:
                print 'key为None，该基因没有对应的pathway'
        print 'pathTable finished'

    def get_pic(self, r_path, map_path, path, kos_path, png_path, html_path, png_bgcolor):
        '''
        画通路图
        '''
        map_id = re.sub('ko', 'map', path)
        pid = re.sub('map', 'ko', path)
        map_html = KeggHtml()
        map_html.color_bg[0] = png_bgcolor
        map_html.run(self.html_path + '/' + map_id + '.html', html_path, path + '.png', self.ko2gene)
        ko_list = [set(self.ko2gene.keys())]
        html_mark = html_path + '.mark'
        map_html.run_html_mark(self.html_path + '/' + map_id + '.html', html_mark, path + '.png', self.ko2gene, ko_list)
        cmd = '{} {} {} {} {} {} {}'.format(
            r_path,
            map_path,
            path,
            kos_path,
            png_path,
            self.html_path + '/' + map_id + '.kgml',
            self.html_path + '/' + map_id + '.png'
        )
        try:
            subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError:
            print '{}画图出错'.format(path)
            os.system('cp {} {}'.format('pathway.png', png_path))

    def keggLayer(self, pathway_table, layerfile):
        '''
        输入pathway_table.xls，获取分类信息文件
        '''
        f = open(pathway_table)
        d = dict()
        ko = dict()
        for record in islice(f, 1, None):
            iterm = record.strip('\n').split('\t')
            ko_list = iterm[5].split(';')
            pid = re.sub('map', 'ko', iterm[0])
            result = self.ko_coll.find_one({'pathway_id': {'$in': [pid]}})
            if result:
                pids = result['pathway_id']
                layer = False
                for index, i in enumerate(pids):
                    if i == pid:
                        category = result['pathway_category'][index]
                        layer_1st = category[0]
                        layer_2nd = category[1]
                        layer = True
                if layer:
                    if ko.has_key(layer_1st):
                        if ko[layer_1st].has_key(layer_2nd):
                            for k in ko_list:
                                ko[layer_1st][layer_2nd].append(k)
                        else:
                            ko[layer_1st][layer_2nd] = ko_list
                    else:
                        ko[layer_1st] = dict()
                        ko[layer_1st][layer_2nd] = ko_list
        with open(layerfile, 'w+') as k:
            for i in ['Metabolism','Genetic Information Processing','Environmental Information Processing','Cellular Processes','Organismal Systems', 'Human Diseases','Drug Development']:
                if ko.has_key(i):
                    for j in ko[i]:
                        ko[i][j] = list(set(ko[i][j]))
                        line = i + '\t' + j + '\t' + str(len(ko[i][j])) + '\t' + ';'.join(ko[i][j]) + '\n'
                        k.write(line)

    def getPic(self, pidpath, pathwaydir, image_magick=None):
        """
        输入文件pid.txt
        输出文件夹pathways
        软件：image_magick：
        功能：作图与将pdf转为png
        目录：~/app/program/ImageMagick/bin/convert
        """
        fs = gridfs.GridFS(self.mongodb)
        f = open(pidpath)
        if not os.path.exists(pathwaydir):
            os.makedirs(pathwaydir)
        for i in f:
            if i:
                i = i.strip('\n').split('\t')
                pid = i[0]
                koid = i[1].split(';')
                l = list()
                kgml_path = os.path.join(os.getcwd(), "pathway.kgml")
                png_path = os.path.join(os.getcwd(), "pathway.png")
                if os.path.exists(kgml_path) and os.path.exists(png_path):
                    os.remove(kgml_path)
                    os.remove(png_path)
                with open("pathway.kgml", "w+") as k, open("pathway.png", "w+") as p:
                    result = self.png_coll.find_one({"pathway_id": pid})
                    if result:
                        kgml_id = result['pathway_ko_kgml']
                        png_id = result['pathway_ko_png']
                        k.write(fs.get(kgml_id).read())
                        p.write(fs.get(png_id).read())
                p_kgml = KGML_parser.read(open("pathway.kgml"))
                p_kgml.image = png_path
                for ortholog in p_kgml.orthologs:
                    for g in ortholog.graphics:
                        g.bgcolor = colors.Color(alpha=0)
                for ko in koid:
                    for degree in p_kgml.entries.values():
                        if re.search(ko, degree.name):
                            l.append(degree.id)
                    for n in l:
                        for graphic in p_kgml.entries[n].graphics:
                            graphic.fgcolor = '#CC0000'
                canvas = KGMLCanvas(p_kgml, import_imagemap=True, label_compounds=True,
                                    label_orthologs=False, label_reaction_entries=False,
                                    label_maps=False, show_maps=False, draw_relations=False, show_orthologs=True,
                                    show_compounds=False, show_genes=False,
                                    show_reaction_entries=False)
                pdf = pathwaydir + '/' + pid + '.pdf'
                png = pathwaydir + '/' + pid + '.png'
                canvas.draw(pdf)
                if image_magick:
                    cmd = '{} -flatten -quality 100 -density 130 -background white {} {}'.format(
                        image_magick, pdf, png
                    )
                    try:
                        subprocess.check_output(cmd, shell=True)
                    except subprocess.CalledProcessError:
                        print '图片格式pdf转png出错'
        print "getPic finished"

    def run(self, r_path, map_path, blast_xml, kegg_ids, kegg_table, pidpath, pathwaydir, pathway_table, layerfile, taxonomy=None, link_bgcolor="green", png_bgcolor="#00CD00", image_magick=None, html_path="/mnt/ilustre/users/sanger-dev/app/database/KEGG/kegg_2017-05-01/kegg/pathway/map", known_ko=None):
        '''
        blast_xml存在对比对到kegg库的xml文件进行kegg注释统计，kegg_ids存在，对客户上传的kegg注释文件进行kegg注释统计
        '''
        if blast_xml:
            print 'INFO: start calling {}'.format(self.pathSearch)
            self.pathSearch(blast_xml, kegg_table, taxonomy)
        if kegg_ids:
            print 'INFO: start calling {}'.format(self.pathSearch_upload)
            self.pathSearch_upload(kegg_ids, kegg_table, taxonomy)
        if known_ko and os.path.isfile(known_ko):
            print 'INFO: start calling {}'.format(self.merge_known)
            self.merge_known(known_ko, kegg_table)
        self.html_path = html_path
        print 'INFO: start calling {}'.format(self.pathTable)
        if png_bgcolor[0] != '#':
            png_bgcolor = '#{}'.format(png_bgcolor)
        self.pathTable(r_path, map_path, kegg_table, pathway_table, pidpath, link_bgcolor, png_bgcolor, pathwaydir, image_magick)
        print 'INFO: start calling {}'.format(self.keggLayer)
        self.keggLayer(pathway_table, layerfile)
        print 'INFO: succeed in running {}'.format(self)

if __name__ == '__main__':
    print 'INFO: start checking incoming parameters before modification'
    for n, i in enumerate(sys.argv):
        print 'sys.argv[{}] = {}'.format(n, i)
    print 'INFO: start modifying sys.argv'
    if sys.argv[3] == 'None':
        sys.argv[3] = None
    if sys.argv[4] == 'None':
        sys.argv[4] = None
    if sys.argv[10] == 'None':
        sys.argv[10] = None
    if sys.argv[11] == 'None':
        sys.argv[11] = 'green'
    if sys.argv[12] == 'None':
        sys.argv[12] = '00CD00'
    if sys.argv[13] == 'None':
        sys.argv[13] = None
    if len(sys.argv) > 15 and sys.argv[15] != 'None':
        known_ko = sys.argv[15]
    else:
        known_ko = None
    print 'INFO: start checking incoming parameters before modification'
    for n, i in enumerate(sys.argv):
        print 'sys.argv[{}] = {}'.format(n, i)
    kegg_anno = KeggAnnotation()
    print 'INFO: start instancing {} as {}'.format(KeggAnnotation, kegg_anno)
    kegg_anno.run(
        r_path=sys.argv[1],
        map_path=sys.argv[2],
        blast_xml=sys.argv[3],
        kegg_ids=sys.argv[4],
        kegg_table=sys.argv[5],
        pidpath=sys.argv[6],
        pathwaydir=sys.argv[7],
        pathway_table=sys.argv[8],
        layerfile=sys.argv[9],
        taxonomy=sys.argv[10],
        link_bgcolor=sys.argv[11],
        png_bgcolor=sys.argv[12],
        image_magick=sys.argv[13],
        html_path=sys.argv[14],
        known_ko=known_ko
    )
