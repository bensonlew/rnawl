# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
#from biocluster.config import config
from collections import namedtuple, defaultdict
import re
import os
import sys
import shutil


class AncestorAgeAgent(Agent):
    def __init__(self, parent):
        super(AncestorAgeAgent, self).__init__(parent)
        options = [
            {"name": "vcf", "type": "infile", "format": "dna_gmap.vcf"},
            {"name":"sample1","type":"infile","format": "align.bwa.bam"},#用作共祖时间分析的bam文件名
            # {"name":"sample2","type":"string"},#用作共祖时间分析的bam文件名2
            {"name": "bam_list", "type": "infile",
                "format": "denovo_rna_v2.bamlist"},
            {"name": "haplogroup", "type": "string"},
            {"name": "out_group", "type": "int", "default": 3},
            {"name": "output", "type": "string", "default": "ytree"},
        ]
        self.add_option(options)

    def check_options(self):
        '''
        参数检查
        '''
        if not self.option("vcf").is_set:
            raise OptionError("请检查vcf文件是否生成")
        if not self.option("bam_list"):
            raise OptionError("必须设置输入bam.list文件")
        if not self.option("haplogroup"):
            raise OptionError("必须输入样本所属单倍型")
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["ageout.txt", "txt", "结果输出文件"],
        ])
        super(AncestorAgeAgent, self).end()


class AncestorAgeTool(Tool):

    def __init__(self, config):
        super(AncestorAgeTool, self).__init__(config)

        #self.ageout = open(os.path.join(self.work_dir, "ageout.txt"), 'w')
        self.s1 = os.path.basename(self.option('sample1').prop["path"]).split('/')[-1].split('.')[0]
        # self.s2 = os.path.basename(self.option('sample2').prop["path"]).split('/')[-1].split('.')[0]
        self.ageout = open(os.path.join(self.work_dir, "{}_ancestor_age.txt".format(self.s1)), 'w')
        self._version = '1.0'
        self.software_dir = self.config.SOFTWARE_DIR
        self.python = 'program/Python/bin/python'
        self.script = os.path.join(
            self.software_dir, 'bioinfo/tool_lab/ancestor_age/ancestor_age.py')
        self.file = {
            'posdb': os.path.join(self.software_dir, 'database/Tool_lab/ysource/posdb'),
            'badpos': os.path.join(self.software_dir, 'database/Tool_lab/ysource/badpos'),
            'tree': os.path.join(self.software_dir, 'database/Tool_lab/ysource/tree'),

        }

    def run(self):
        '''
        运行
        '''
        super(AncestorAgeTool, self).run()
        self.run_Ancestor_age()
        self.set_output()
        self.end()

    def run_Ancestor_age(self):
        '''
        python ancestor_age.py -t -p -b -v -l -g -d -o -u
        '''
        posdb = self.file['posdb']
        tree = self.file['tree']
        badpos = self.file['badpos']
        vcf = self.option('vcf').prop['path']
        bam_list = self.option('bam_list').prop['path']
        haplogroup = self.option('haplogroup')
        outgroup = self.option('out_group')
        cmd = '{} {}'.format(self.python, self.script)
        cmd += ' -t {} -p {} -b {}'.format(tree, posdb, badpos)
        cmd += ' -v {} -l {} -g {} -o {}'.format(
            vcf, bam_list, haplogroup, outgroup)
        cmd += ' -u {}'.format(os.path.join(self.work_dir, "{}_ancestor_age.txt".format(self.s1)))
        self.logger.info(cmd)
        self.logger.info("开始计算共祖时间")
        command1 = self.add_command("ancestor_age", cmd).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("运行ancestor_age完成")
        else:
            self.logger.info("运行错误，请重新检查输入参数")

    def set_output(self):
        '''
        将结果文件赋值到output文件夹下面
        :return:
        '''
        if len(self.output_dir) > 0:
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        self.logger.info("设置结果目录")
        try:
            os.link(os.path.join(self.work_dir, "{}_ancestor_age.txt".format(self.s1)),
                    os.path.join(self.output_dir, "{}_ancestor_age.txt".format(self.s1)))
        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))

    # def linkfile(self):
    #     age_file = os.path.join(self.work_dir, 'ageout.txt')
    #     link = self.output_dir + '/ageout.txt'
    #     if os.path.exist(link):
    #         os.remove(link)
    #     os.link(age_file, link)

    # def set_db(self):
    #     self.logger.info("开始导表")
    #     api_ancestorage = self.api.api("tool_lab.ancestor_age")
    #     api_ancestorage.add_ancestorage_detail(self.option('main_id'), self.output_dir)
    #     self.logger.info("导表结束")


"""
class AncestorAgeTool(Tool):
    
    def __init__(self, config):
        super(AncestorAgeTool, self).__init__(config)
        
        self.ageout = open(os.path.join(self.work_dir, "ageout.txt"),'w')
        self._version = '1.0'
        self.software_dir = self.config.SOFTWARE_DIR
        self.file = {
            'posdb': os.path.join(self.software_dir, 'database/Tool_lab/ysource/posdb'),
            'badpos': os.path.join(self.software_dir, 'database/Tool_lab/ysource/badpos'),
            'tree': os.path.join(self.software_dir, 'database/Tool_lab/ysource/tree'),
            
        }

        
    def run(self):
        '''
        运行
        '''
        super(AncestorAgeTool, self).run()
        self.run_Ancestor_age()

    def run_Ancestor_age(self):
        Clade = namedtuple('Clade', ['id', 'name', 'level', 'parent', 'children', 'muts', 'poss', 'alts','ads','old','time','hpd','envidence','population'])
        self.ageout.write('Begin>>')
        self.logger.info("Begin>>")
        posdb = self.file['posdb']
        tree = self.file['tree']
        badpos = self.file['badpos']
        tr = Tree(tree ,posdb ,badpos ,self.option('vcf'),self.option('bam_list'),self.option('haplogroup'),self.option('output'),self.option('dep'),self.option('outgroup'))
        self.logger.info("{}".format(tr))
        self.ageout.write("\n{}".format(tr))
        self.ageout.close()

    def set_output(self):
        '''
        将结果文件赋值到output文件夹下面
        :return:
        '''
        self.logger.info("设置结果目录")
        try:
            os.link(os.path.join(self.work_dir, "ageout.txt"), os.path.join(self.output_dir, "ageout.txt"))
        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))        



    def linkfile(self):
        age_file = os.path.join(self.work_dir, 'ageout.txt') 
        link = self.output_dir + '/ageout.txt'
        if os.path.exist(link):
            os.remove(link)
        os.link(age_file, link)

    def set_db(self):
        self.logger.info("开始导表")
        api_scatter = self.api.api("tool_lab.ancestor_age_api")
        api_scatter.add_scater_detail(self.option('main_id'), self.output_dir)
        self.logger.info("导表结束")


class PosDb(object):
    '''posdb
    '''
    
    field = {'name','pos','ref','alt',}
    def __init__(self, db_file, bad_file=None):
        self._pos_dic = {} # pos: Mut()
        self.ageout = open(os.path.join(self.work_dir, "ageout.txt"),'a')
        self._pos_class = namedtuple('Pos', ['pos','name','ref','alt'])
        self.name_pos = {}  # name: pos
        self.badpos = {}
        self.onref = {}
        self.pos_clade = {}

        if bad_file:
            self._read_badpos(bad_file)
        self._read_dbfile(db_file)
        self.poss = self._pos_dic.keys()

    def _read_dbfile(self, db_file):
        self.ageout.write("\nstart read db...")
        self.logger.info("start read db...")
        with open(db_file, "r") as db:
            head = db.readline()
            for line in db:
                line = re.sub('"|\'','',line)
                tabs = line.rstrip('\n').split('\t')

                p, ns, r, a, clade, onref = tabs[:6]

                if r in ['']:
                    r = "?"
                if a == '':

                    a = '?'

                if len(tabs)>6:
                    unreliable = tabs[6]
                    if unreliable in ['b','t']:
                        self.badpos[p] = unreliable
                self.pos_clade[p] = clade

                self._pos_dic[p] = self._pos_class(p, ns, r, a)

                for n in ns.split('/'):

                    self.name_pos[n] = p
            self.ageout.write("\n...ok")
            self.logger.info("...ok")

    def _read_badpos(self, bad_file):
        self.ageout.write("\nstart read badpos...")
        self.logger.info("start read badpos...")
        with open(bad_file, "r") as db:
            for line in db:
                tabs = line.rstrip('\n').split('\t')
                self.badpos[tabs[0]]=tabs[1]

        self.ageout.write("\n...ok")
        self.logger.info("...ok")

    def add_badpos(self, pos, tag):
        if pos in self.badpos.keys():
            pass
        else:
            self.badpos[pos] = tag

    def add_onref(self,pos,tag):
        if pos in self.onref.keys():
            if self.onref[pos] != tag:
                self.logger.info(">>>>onref conflict:{}".format(pos))
            else:
                self.logger.info((">>>onref recurrent:{}".format(pos)))
        else:
            self.onref[pos] = tag

    def __getitem__(self, pos):
        return self._pos_dic[pos]

    def update_field(self, pos, field, value):

        repd = {field:value}
        self._pos_dic[pos] = self._pos_dic[pos]._replace(**repd)

    def update_fields(self, pos, fields, values):

        repd = dict(zip(fields,values))
        self._pos_dic[pos] = self._pos_dic[pos]._replace(**repd)

    def __repr__(self):
        return('This is a posdb: {} snp!'.format(len(self.poss)))

    def add_pos(self, value):
        if value[0] not in self._pos_dic.keys():
            self._pos_dic[value[0]] = self._pos_class(*value)


    def add_pos_clade(self, pos, clade):
        if pos not in self.pos_clade.keys():
            self.pos_clade[pos] = clade
        if self.pos_clade[pos] !='':
            if self.pos_clade[pos] != clade:
                self.pos_clade[pos] = self.pos_clade[pos] + "," + clade
        else:
            self.pos_clade[pos] = clade

class Vcf(object):
    '''vcf
    '''
    def __init__(self, vcf_file, posdb, bam_list = None, filter_bad=False, dep = 5):
        self.ageout = open(os.path.join(self.work_dir, "ageout.txt"),'a')
        self.posdb = posdb
        self.clades = {}
        self.pos_names = []
        self.filter_bad = filter_bad
        self.samples_name = []
        self.dep = dep
        self._read_vcf(vcf_file,bam_list)

    def _read_vcf(self, vcf_file, bam_list):
        #读取多样本mpileup后的vcf文件，获取样本的突变位点list
        self.ageout.write("\nstart read vcf...")
        self.logger.info("start read vcf...")

        if bam_list:
            with open(bam_list, 'r') as bl:
                for l in bl:
                    s = l.split('/')[-1].split('.')[0]
                    self.samples_name.append(s)
        with open(vcf_file, 'r') as vcf:
            while 1:
                line = vcf.readline()
                if line[0:2] == '##':
                    pass
                elif line[0] == '#':
                    if bam_list:
                        pass
                    else:
                        self.samples_name = line.rstrip().split('\t')[9:]
                    for sn in self.samples_name:

                        self.clades[sn] = Clade(-1, sn, 0, 'unknow', [], [], [],[],[],'','','','','')
                    break

            while 1:

                line = vcf.readline()
                if not line:
                    break
                fields = line.rstrip().split('\t')
                if fields[7].split(";")[0] == 'INDEL':
                    continue
                pos, ref, alts, formats, samples_gt = fields[1], fields[3], fields[4], fields[8].split(':'), fields[9:]

                gt_i, ad_i =  formats.index('GT'), formats.index('AD')

                ref_alts = alts.split(',')
                ref_alts.insert(0,ref)
                s_gtypes = []

                for sgt in samples_gt:

                    gt, ad = sgt.split(':')[gt_i], sgt.split(':')[ad_i]

                    gtype = self._check_genotype( ref_alts,gt, ad.split(','))

                    s_gtypes.append(gtype)

                if self._pos_filter(pos, s_gtypes):  # 粗筛过滤

                    onref = self.posdb.onref[pos] if pos in self.posdb.onref.keys() else ''
                    self.logger.info("process:"+ pos)
                    for sgt in samples_gt:
                        si = samples_gt.index(sgt)

                        flag = 0

                        if s_gtypes[si][0] in ['1','2','3'] and onref != "*": # 同alt
                            flag = 1
                        else:
                            if s_gtypes[si][0] == '0' and onref == "*":
                                flag = 1

                        if flag == 1:

                            self.clades[self.samples_name[si]].poss.append(pos)

                    self.pos_names.append(pos)


        self.ageout.write("\n...ok")
        self.logger.info("...ok")

    def _check_genotype(self, alts, gt, ad):
        #检查突变类型

        alt_i = 0
        ads = [int(x) for x in ad]
        alt_i = ads.index(max(ads))

        gtype = {
            './.': '-',
            '0/0': '0'+str(ads[alt_i]),
            '1/1': '1'+str(ads[alt_i]),
            '2/2': '2'+str(ads[alt_i]),
            '3/3': '3'+str(ads[alt_i]),
            '0/1': '?',
            '0/2': '?',
            '0/3': '?',
            '1/2': str(alt_i)+str(ads[alt_i]),
            '1/3': str(alt_i)+str(ads[alt_i]),
            '2/3': str(alt_i)+str(ads[alt_i]),
        }
        return gtype[gt]


    def _pos_filter(self, pos, gtypes):
        # 变异位点粗筛, 计算年份时不使用
        if pos in self.posdb.badpos.keys():
            return False
        if gtypes.count('-') > 0:
            return False
        if gtypes.count('?') > 0:
            return False
        def filter_cov(x):
            flag = 1 if x >= self.dep else 0
            return flag
        flags = [filter_cov(int(x[1:])) for x in gtypes]
        if sum(flags) < len(flags):
            return False
        return True

    def __repr__(self):
        return('This is a Vcf: {} samples, {} pos !'.format(len(self.clades), len(self.pos_names)))

class Tree(object):
    '''tree
    '''

    def __init__(self, tree_file, posdb_file=None, badpos_file=None,vcf_file = None,bam_list=None,haplogroup = None, output = 'ytree', dep = 5,outgroup =3):
        self.ageout = open(os.path.join(self.work_dir, "ageout.txt"),'a')
        self.output = output
        self.posdb = PosDb(posdb_file,badpos_file)
        ageout.write("\n%s"%self.posdb)
        self.logger.info(self.posdb)
        self._clades = {}  # name:clade
        self.pos_names = []
        self.tree_file = tree_file
        self._read_tree()
        self.vcf = None
        if vcf_file !=None:
            self.vcf = Vcf(vcf_file, self.posdb, bam_list,dep=dep)  # , filter_bad=calc_age
            self.ageout.write("\n%s"%self.vcf)
            self.logger.info(self.vcf)
            self._clades.update(self.vcf.clades)
            self.pos_names = list(set(self.pos_names + self.vcf.pos_names))
        self.calc_time(haplogroup,int(outgroup))


    def _read_tree(self):
        # 按行读取树文件，获取节点的突变类型
        self.ageout.write('\nstart read tree...')
        self.logger.info("start read tree...")
        with open(self.tree_file, "r") as tree:
            tree.readline() # title行
            root_line = tree.readline().rstrip()  # root行
            l = 0
            time,hpd,envidence,population = root_line.split('\t')[0:4]
            self._clades['0root'] = Clade(l,"0root", 0, None, [], [], [],[],[],'',time,hpd,envidence,population)
            last_clade = self._clades['0root']
            for line in tree:
                l += 1
                time,hpd,envidence,population = line.rstrip().split('\t')[0:4]
                tr = '\t'.join(line.rstrip().split('\t')[4:])
                if not tr:
                    continue
                m = re.match(r'^(\s*)([^-\s]+)-(.*)', tr)
                if not m:
                    self.ageout.write("\nline {} format not match!{}".format(l, line))
                    self.logger.info("line {} format not match{}".format(1, line))
                    continue
                level, name, muts = (len(m.group(1)), m.group(2), re.split(r',\s*', m.group(3)))
                while level <= last_clade.level:
                    if last_clade.parent is None:
                        break
                    else:
                        last_clade = self._clades[last_clade.parent]
                poss, alts = [None,None]
                if self.posdb:
                    poss, alts = self._get_pos_by_name(muts,name)
                    self.pos_names += poss
                    alts.extend(last_clade.alts)
                    poss.extend(last_clade.poss)
                    muts.extend(last_clade.muts)

                self._clades[name]  = Clade(l, name, level, last_clade.name, [], muts, poss, alts,[],'',time,hpd,envidence,population)
                self._clades[last_clade.name].children.append(name)
                last_clade = self._clades[name]
        self.ageout.write("\n...ok")
        self.logger.info("...ok")

    def calc_time(self, haplogroup, outgroup):
        self.ageout.write('\n>>>start calc time with haplogroup:'+haplogroup)
        self.logger.info(">>>start calc time with haplogroup:"+haplogroup)
        know,t0 = self._find_closest_time(haplogroup)
        self.ageout.write("\n>>>t0 time:"+know+":"+str(t0))
        self.logger.info(">>>t0 time:"+know+":"+str(t0))
        samples = self.vcf.samples_name
        self.ageout.write("\n>>>outgroups:"+','.join(samples[outgroup-1:]))
        self.logger.info(">>>outgroups:"+','.join(samples[outgroup-1:]))
        outgroup_set = set(self._clades[samples[outgroup-1]].poss).union(*[self._clades[x].poss for x in samples[outgroup:]])
        self.ageout.write("\n>>>samples:"+','.join(samples[:outgroup-1]))
        self.logger.info(">>>samples:"+','.join(samples[:outgroup-1]))
        sr = samples[0]
        sr_set = set(self._clades[sr].poss)-outgroup_set
        all_anc_set = sr_set.intersection(*[self._clades[sq].poss for sq in samples[1:outgroup-1]])
        all_s_len = len(sr_set-all_anc_set)
        self.ageout.write("\n>>>"+sr+":after filter:"+str(len(sr_set)))
        self.logger.info(">>>"+sr+":after filter:"+str(len(sr_set)))
        for sq in samples[1:outgroup-1]:
            sq_set = set(self._clades[sq].poss)-outgroup_set
            self.ageout.write("\n>>>"+sq+":after filter:"+str(len(sq_set)))
            self.logger.info(">>>"+sq+":after filter:"+str(len(sq_set)))
            anc_set = sr_set.intersection(sq_set)
            anc_len = len(anc_set)
            self.ageout.write("\n>>>"+"anc_len:"+str(anc_len))
            self.logger.info(">>>"+"anc_len:"+str(anc_len))
            sr_len = len(sr_set-anc_set)
            self.ageout.write("\n>>>"+"sr_len:"+str(sr_len))
            self.logger.info(">>>"+"sr_len:"+str(sr_len))
            sq_len = len(sq_set-anc_set)
            all_s_len += len(sq_set-all_anc_set)
            self.ageout.write("\n>>>"+"sq_len:"+str(sq_len))
            self.logger.info(">>>"+"sq_len:"+str(sq_len))
            s_len = 0.5*(sr_len+sq_len)
            anc_time = t0*s_len/(s_len+anc_len)
            self.ageout.write("\n>>>"+sr+'-'+sq +':'+ str(anc_time) + '\n')
            self.logger.info(">>>"+sr+'-'+sq +':'+ str(anc_time) + '\n')
        all_s_len = all_s_len/(outgroup-1)
        all_anc_len = len(all_anc_set)
        all_anc_time = t0*all_s_len/(all_s_len+all_anc_len)
        self.ageout.write("\n>>>All samples:"+str(all_anc_time) +'\n')
        self.logger.info(">>>All samples:"+str(all_anc_time) +'\n')

    def _pos_filter(self, clade, dep=2):
        # 变异位点粗筛
        del_set = []
        for pos in self._clades[clade].poss:
            flag = 0
            if len(self._clades[clade].ads) and self._clades[clade].ads[self._clades[clade].poss.index(pos)] < dep:
                self.ageout.write("\nfilter"+str(dep)+":"+pos+":"+str(self._clades[clade].ads[self._clades[clade].poss.index(pos)]))
                self.logger.info("filter"+str(dep)+":"+pos+":"+str(self._clades[clade].ads[self._clades[clade].poss.index(pos)]))
                flag =1
            else:
                if pos in self.posdb.badpos.keys():
                    flag = 1

            if flag:
                del_set.append(pos)

        return set(self._clades[clade].poss)-set(del_set)

    def _find_closest_time(self,clade):

        if self._clades[clade].time =='':
            parent = self._clades[clade].parent
            return self._find_closest_time(parent)

        return [clade,float(self._clades[clade].time)]

    def __repr__(self):
        return('This is a Tree: {} clades, {} mutations'.format(len(self._clades), len(self.pos_names)))

    def __iter__(self):
        for clade in self._clades:
            yield clade
        return

    def __getitem__(self, name):
        return self._clades[name]

    def _get_pos_by_name(self, names,haplogroup):

        poss = []
        alts = []

        for name in names:

            pos = "none"
            alt = "?"
            nms = name.split('/')
            if name[-1] == "#":
                continue
            elif name[-1] == '&':
                nm  = nms[0] if len(nms)>1 else nms[0][:-1]
                if nm in self.posdb.name_pos.keys():
                    pos = self.posdb.name_pos[nm]
                    pd = self.posdb[pos]
                    alt = pd.ref

                    self.posdb.add_onref(pos,'*')
                    self.posdb.add_pos_clade(pos,haplogroup)

                else:
                    pos = nm
                    alt = '?'

                    self.posdb.add_onref(pos,'*')
                    self.posdb.add_pos([pos,name[:-1],'?','?'])
                    self.posdb.add_pos_clade(pos,haplogroup)
                    self.ageout.write("\n"+nm+" not found!")
                    self.logger.info(nm+" not found!")

            elif name[-2] == '=':
                nm  = nms[0] if len(nms)>1 else nms[0][:-2]
                if nm in self.posdb.name_pos.keys():
                    pos = self.posdb.name_pos[nm]
                    self.posdb.add_pos_clade(pos,haplogroup)
                    self.posdb.add_onref(pos,'')
                    alt = name[-1]

                else:
                    pos = nm
                    alt = name[-1]

                    self.posdb.add_pos([pos,name[:-2],'N',alt])
                    self.posdb.add_pos_clade(pos,haplogroup)
                    self.posdb.add_onref(pos,'')
                    self.ageout.write("\n"+nm+" not found!")
                    self.logger.info(nm+" not found!")
            else:
                nm = nms[0]
                if nm in self.posdb.name_pos.keys():
                    pos = self.posdb.name_pos[nm]
                    self.posdb.add_pos_clade(pos,haplogroup)
                    self.posdb.add_onref(pos,'')
                    pd = self.posdb[pos]

                    if self.posdb.onref[pos]=='*':
                        alt = pd.ref[0]
                    else:
                        alt = pd.alt[0]

                else:
                    pos =nm
                    alt = '?'

                    self.posdb.add_pos([pos,name,'N','?'])
                    self.posdb.add_pos_clade(pos,haplogroup)
                    self.posdb.add_onref(pos,'')
                    self.ageout.write("\n"+nm+" not found!")
                    self.logger.info(nm+" not found!")

            poss.append(pos)
            alts.append(alt)
        return [poss,alts]
"""
