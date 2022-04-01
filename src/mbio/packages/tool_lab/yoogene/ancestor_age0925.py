# -*- coding: utf-8 -*-
# __author__ = 'yuguo'

import re, os, sys
import argparse
# import subprocess
# import pickle
from collections import namedtuple, defaultdict
# from ete3 import PhyloTree,TreeStyle, TextFace

Clade = namedtuple('Clade', ['id', 'name', 'level', 'parent', 'children', 'muts', 'poss', 'alts','ads','old','time','hpd','envidence','population'])

class PosDb(object):
    '''posdb
    '''
    # field = {'name','pos','ref','alt','onref','unreliable','haplogroup'}
    field = {'name','pos','ref','alt',}
    def __init__(self, db_file, bad_file=None):
        self._pos_dic = {} # pos: Mut()
        # self.fields = []
        # self._pos_class = namedtuple('Pos', ['pos','name','ref','alt','onref','unreliable','haplogroup'])
        # self._pos_class = namedtuple('Pos', ['pos','name','ref','alt','onref'])
        self._pos_class = namedtuple('Pos', ['pos','name','ref','alt'])
        self.name_pos = {}  # name: pos
        self.badpos = {}
        self.onref = {}
        self.pos_clade = {}
        # self.name_onref = {} # pos: name+onref
        if bad_file:
            self._read_badpos(bad_file)
        self._read_dbfile(db_file)
        self.poss = self._pos_dic.keys()

    def _read_dbfile(self, db_file):
        print("start read db...")
        agetext.write("\nstart read db...")
        with open(db_file, "r") as db:
            head = db.readline()
            # self.fields = db.readline().split('\t')
            # self._pos_class = namedtuple('Pos', db.readline().rstrip().split('\t'))
            # if not self.field <=set(self._pos_class._fields):
            #     raise Exception("Pos db title fields not match: {}".format(str(self.field - self._pos_class._fields)))
            for line in db:
                line = re.sub('"|\'','',line)
                tabs = line.rstrip('\n').split('\t')
                # print(tabs)
                # onref = tabs[4] if len(tabs)>4 else ''
                # if len(tabs)<4:
                #     tabs.append('?')
                p, ns, r, a, clade, onref = tabs[:6]
                # if r in ['ins', 'del'] or a in ['ins', 'del']:
                #     continue
                if r in ['']:
                    r = "?"
                if a == '':
                    # a = "-" if onref == 'indel' else '?'
                    a = '?'
                # self.onref[p] = ''
                if len(tabs)>6:
                    unreliable = tabs[6]
                    if unreliable in ['b','t']:
                        self.badpos[p] = unreliable
                self.pos_clade[p] = clade
                # if unreliable == '':
                #     unreliable = self.badpos[p] if p in self.badpos.keys() else ''
                # else:
                #     self.badpos[p] = unreliable if p not in self.badpos.keys() else self.badpos[p]
                # r = r[0]
                # a = a[0]
                # if len(r) >1 or len(a) >1:
                #     # print(r+'-'+a)
                #     r = r[0]
                #     a = a[0]
                # self._pos_dic[p] = self._pos_class(p, ns, r.upper(), a.upper())
                # self._pos_dic[p] = self._pos_class(p, ns, r, a, onref,unreliable,[])
                self._pos_dic[p] = self._pos_class(p, ns, r, a)
                # self.poss[p] = ns
                for n in ns.split('/'):
                    # pc = self._pos_class(p, n, r.upper(), a.upper())
                    # if n in self._pos_dic.keys():
                    #     print('WARNING: snp {} present in two position, please check!'.format(n))
                    self.name_pos[n] = p
            print("...ok")
            agetext.write("\n...ok")

    def _read_badpos(self, bad_file):
        print("start read badpos...")
        agetext.write("\nstart read badpos...")
        with open(bad_file, "r") as db:
            for line in db:
                tabs = line.rstrip('\n').split('\t')
                self.badpos[tabs[0]]=tabs[1]
                # if tabs[0] in self._pos_dic.keys():
                #     self.update_field(tabs[0],'unreliable',tabs[1])
                # else:
                #     self._pos_dic[tabs[0]] = self._pos_class(tabs[0], '', '', '', '',tabs[1],[])
        print("...ok")
        agetext.write("\n...ok")

    def add_badpos(self, pos, tag):
        if pos in self.badpos.keys():
            pass
        else:
            self.badpos[pos] = tag

    def add_onref(self,pos,tag):
        if pos in self.onref.keys():
            if self.onref[pos] != tag:
                print(">>>>onref conflict:{}".format(pos))
        else:
            self.onref[pos] = tag

    def __getitem__(self, pos):
        return self._pos_dic[pos]

    def update_field(self, pos, field, value):
        # lpos = list(self._pos_dic[pos])
        # lpos[self._pos_class._fields.index(field)] = value
        # self._pos_dic[pos] = self._pos_class(*lpos)
        repd = {field:value}
        self._pos_dic[pos] = self._pos_dic[pos]._replace(**repd)

    def update_fields(self, pos, fields, values):
        # cla = list(self._clades[key])
        repd = dict(zip(fields,values))
        self._pos_dic[pos] = self._pos_dic[pos]._replace(**repd)

    def __repr__(self):
        return('This is a posdb: {} snp!'.format(len(self.poss)))

    def add_pos(self, value):
        if value[0] not in self._pos_dic.keys():
            self._pos_dic[value[0]] = self._pos_class(*value)
            # self.onref[value[0]] = onref
            # self.poss.append(value[0])

    def add_pos_clade(self, pos, clade):
        if pos not in self.pos_clade.keys():
            self.pos_clade[pos] = clade
        if self.pos_clade[pos] !='':
            if self.pos_clade[pos] != clade:
                self.pos_clade[pos] = self.pos_clade[pos] + "," + clade
        else:
            self.pos_clade[pos] = clade

    def write_posdb(self):
        '''生成基于树更新clade的posdb和bed文件'''
        pbed = open("posdb_update.bed", 'w')
        pdb = open("posdb_update.xls", 'w')
            # pdb.write('\xef\xbb\xbf'+'\t'.join(['pos','name','ref','alt','onref','haplogroup'])+"\n")
        pdb.write('\t'.join(['pos','name','ref','alt','onref','clade'])+"\n")
        for pos in self.poss:
            onref = self.onref[pos] if pos in self.onref.keys() else ''
            pdb.write('\t'.join(list(self._pos_dic[pos]))+'\t'+self.pos_clade[pos]+'\t'+onref+'\n')
            pbed.write('\t'.join(['chrY',str(int(pos)-1),str(pos)])+"\n")
        pbed.close()
        pdb.close()
        


class Vcf(object):
    '''vcf
    '''
    def __init__(self, vcf_file, posdb, bam_list = None, filter_bad=False, dep = 5):
        self.posdb = posdb
        self.clades = {}
        self.pos_names = []
        self.filter_bad = filter_bad
        self.samples_name = []
        self.dep = dep
        self._read_vcf(vcf_file,bam_list)
        # self.samples_pos = defaultdict(list)  # sample: [pos]
    def _read_vcf(self, vcf_file, bam_list):
        """读取多样本mpileup后的vcf文件，获取样本的突变位点list"""
        print("start read vcf...")
        agetext.write("\n\nstart read vcf...")
        
        # samples_name = None
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
                        # Clade = namedtuple('Clade', ['id', 'name', 'level', 'parent', 'children', 'muts', 'poss', 'alts'])
                        self.clades[sn] = Clade(-1, sn, 0, 'unknow', [], [], [],[],[],'','','','','')
                    break

            while 1:
            # for line in vcf:
            #     if line[0:2] == '##':
            #         pass
            #     elif line[0] == '#':
            #         samples_name = line.strip().split('\t')[9:]
            #         for sn in samples_name:
            #             # Clade = namedtuple('Clade', ['id', 'name', 'level', 'parent', 'children', 'muts', 'poss', 'alts'])
            #             self.clades[sn] = Clade(0, sn, -1, 'unknow', [], [], [],[])
            #     else:
                # print(line)
                line = vcf.readline()
                if not line:
                    break
                fields = line.rstrip().split('\t')
                if fields[7].split(";")[0] == 'INDEL':
                    continue
                pos, ref, alts, formats, samples_gt = fields[1], fields[3], fields[4], fields[8].split(':'), fields[9:]
                # print(fields)
                # db_alt = self.posdb[pos].alt
                # if not set('GT','AD') <= set(formats.split(':')):
                #     raise Exception("vcf field format not match both GT/AD")
                # else:
                gt_i, ad_i =  formats.index('GT'), formats.index('AD')
                # sample_num = len(samples_gt)
                # pp = 0
                ref_alts = alts.split(',')
                ref_alts.insert(0,ref)
                s_gtypes = []
                # s_ads = []
                for sgt in samples_gt:
                    # si = samples_gt.index(sgt)
                    # print(samples_name[si])
                    gt, ad = sgt.split(':')[gt_i], sgt.split(':')[ad_i]
                    # gtype, vcfsimp, ad = self._check_genotype( ref_alts,gt, ad.split(','))
                    gtype = self._check_genotype( ref_alts,gt, ad.split(','))
                    # vcfsimps.append(vcfsimp)
                    s_gtypes.append(gtype)
                    # s_ads.append(ad)
                if self._pos_filter(pos, s_gtypes):  # 粗筛过滤
                        # pp = 1
                    # mut = self.posdb[pos].name if pos in self.posdb.poss else pos
                    onref = self.posdb.onref[pos] if pos in self.posdb.onref.keys() else ''
                    #print("process:"+ pos)
                    #agetext.write("\n\nprocess:{}".format(pos))
                    for sgt in samples_gt:
                        si = samples_gt.index(sgt)
                        # self.clades[samples_name[si]].muts.append(mut + ':'+ ref+','+alts)
                        flag = 0
                        # if self.filter_bad and s_ads[si] < 2:
                        # if  s_ads[si] < 2:
                        #     flag = 0
                        # else:
                        if s_gtypes[si][0] in ['1','2','3'] and onref != "*": # 同alt
                            flag = 1
                        else:
                            if s_gtypes[si][0] == '0' and onref == "*":
                                flag = 1
                        # if s_gtypes[si] in ['-','?']:
                        #     flag = 0
                        # if self.filter_bad and s_ads[si] < 5:
                        #         flag = 0
                        if flag == 1:
                            # self.clades[self.samples_name[si]].muts.append(mut)
                            self.clades[self.samples_name[si]].poss.append(pos)
                            # self.clades[self.samples_name[si]].ads.append(s_ads[si])
                            # self.clades[self.samples_name[si]].alts.append(ref_alts[int(s_gtypes[si][0])])
                            # print("flag1:"+self.samples_name[si]+":"+ pos)
                            # self.clades[samples_name[si]].muts.append(self.posdb[pos].name)
                    self.pos_names.append(pos)
                    # self.posdb.add_pos([mut, pos, ref, alts,'','',''])
                    # self.samples_pos[samples_name[si]].append(str(pos)+':'+alt)
                # if pp: self.pos_names.append(pos)

        print("...ok")
        agetext.write("\n...ok")
        

    def _check_genotype(self, alts, gt, ad):
        """检查突变类型"""
        # alts.insert(0,ref)
        alt_i = 0
        ads = [int(x) for x in ad]
        alt_i = ads.index(max(ads))
        # cov = '2' if ads[alt_i] > 1 else '1'
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
        # vcfsimp = {
        #     '-':'-', # 0x覆盖
        #     '?':'?', # 杂合
        #     # '1?':'?', # 杂合
        #     # '2?':'?', # 杂合
        #     '01':'o', # 1x覆盖且同ref
        #     '02':'0', # >=2x覆盖且同ref
        #     '11':'|',  # 1x覆盖且同alt
        #     '12':'1',  # >=2x覆盖且同alt
        #     '21':alts[alt_i].lower(),  # 1x覆盖且不同ref或alt
        #     '22':alts[alt_i].upper(),  # >=2x覆盖且不同ref或alt
        # }
        # return [gtype[gt], vcfsimp[gtype[gt]],ads[alt_i]]

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
        self.output = output
        self.posdb = PosDb(posdb_file,badpos_file)
        print(self.posdb)
        agetext.write("\n{}".format(self.posdb))
        self._clades = {}  # name:clade
        self.pos_names = []
        self.tree_file = tree_file
        self._read_tree()
        self.vcf = None
        if vcf_file !=None:
            self.vcf = Vcf(vcf_file, self.posdb, bam_list,dep=dep)  # , filter_bad=calc_age
            print(self.vcf)
            agetext.write("\n{}".format(self.vcf))
            self._clades.update(self.vcf.clades)
            self.pos_names = list(set(self.pos_names + self.vcf.pos_names))
        self.calc_time(haplogroup,int(outgroup))


    def _read_tree(self):
        # 按行读取树文件，获取节点的突变类型
        print('start read tree...')
        agetext.write('start read tree...')
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
                    print("line {} format not match!{}".format(l, line))
                    agetext.write("\nline {} format not match!{}".format(l, line))
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
        print("...ok")
        agetext.write("\n...ok")

    def calc_time(self, haplogroup, outgroup):
        print('>>>start calc time with haplogroup:'+haplogroup)
        agetext.write('>>>start calc time with haplogroup:'+haplogroup)
        know,t0 = self._find_closest_time(haplogroup)
        print(">>>t0 time:"+know+":"+str(t0))
        agetext.write("\n>>>t0 time:"+know+":"+str(t0))
        samples = self.vcf.samples_name
        print(">>>outgroups:"+','.join(samples[outgroup-1:]))
        agetext.write("\n>>>outgroups:"+','.join(samples[outgroup-1:]))
        outgroup_set = set(self._clades[samples[outgroup-1]].poss).union(*[self._clades[x].poss for x in samples[outgroup:]])
        print(">>>samples:"+','.join(samples[:outgroup-1]))
        agetext.write("\n>>>samples:"+','.join(samples[:outgroup-1]))
        sr = samples[0]
        sr_set = set(self._clades[sr].poss)-outgroup_set
        all_anc_set = sr_set.intersection(*[self._clades[sq].poss for sq in samples[1:outgroup-1]])
        all_s_len = len(sr_set-all_anc_set)
        print(">>>"+sr+":after filter:"+str(len(sr_set)))
        agetext.write("\n>>>"+sr+":after filter:"+str(len(sr_set)))
        for sq in samples[1:outgroup-1]:
            sq_set = set(self._clades[sq].poss)-outgroup_set
            print(">>>"+sq+":after filter:"+str(len(sq_set)))
            agetext.write("\n>>>"+sq+":after filter:"+str(len(sq_set)))
            anc_set = sr_set.intersection(sq_set)
            anc_len = len(anc_set)
            print(">>>"+"anc_len:"+str(anc_len))
            agetext.write("\n>>>"+"anc_len:"+str(anc_len))
            sr_len = len(sr_set-anc_set)
            print(">>>"+"sr_len:"+str(sr_len))
            agetext.write("\n>>>"+"sr_len:"+str(sr_len))
            sq_len = len(sq_set-anc_set)
            all_s_len += len(sq_set-all_anc_set)
            print(">>>"+"sq_len:"+str(sq_len))
            agetext.write("\n>>>"+"sq_len:"+str(sq_len))
            s_len = 0.5*(sr_len+sq_len)
            anc_time = t0*s_len/(s_len+anc_len)
            print(">>>"+sr+'-'+sq +':'+ str(anc_time) + '\n')
            agetext.write("\n>>>"+sr+'-'+sq +':'+ str(anc_time) + '\n')
            anc_pos = '\t'.join(anc_set)
            print(">>>anc_pos:" + anc_pos + '\n')
            agetext.write("\n>>>anc_pos:" + anc_pos + '\n')
            sr_pos = '\t'.join(sr_set - anc_set)
            print(">>>sr_pos:" + sr_pos + '\n')
            agetext.write("\n>>>sr_pos:" + sr_pos + '\n')
            sq_pos = '\t'.join(sq_set - anc_set)
            print(">>>sq_pos:" + sq_pos + '\n')
            agetext.write("\n>>>sq_pos:" + sq_pos + '\n')
        all_s_len = all_s_len/(outgroup-1)
        all_anc_len = len(all_anc_set)
        all_anc_time = t0*all_s_len/(all_s_len+all_anc_len)
        print(">>>All samples:"+str(all_anc_time) +'\n')
        agetext.write("\n>>>All samples:"+str(all_anc_time) +'\n')

    def _pos_filter(self, clade, dep=2):
        # 变异位点粗筛
        del_set = []
        for pos in self._clades[clade].poss:
            flag = 0
            if len(self._clades[clade].ads) and self._clades[clade].ads[self._clades[clade].poss.index(pos)] < dep:
                print("filter"+str(dep)+":"+pos+":"+str(self._clades[clade].ads[self._clades[clade].poss.index(pos)]))
                agetext.write("\nfilter"+str(dep)+":"+pos+":"+str(self._clades[clade].ads[self._clades[clade].poss.index(pos)]))
                flag =1
            else:
                if pos in self.posdb.badpos.keys():
                    flag = 1
                    # unreliable = self.posdb[pos].unreliable  # 检查位点库中位点的unreliable属性
                    # if unreliable in ['b','t']:
                    #     flag =1
            if flag:
                del_set.append(pos)
                # if len(self._clades[clade].ads):
                #     del self._clades[clade].ads[self._clades[clade].poss.index(pos)]
                # del self._clades[clade].alts[self._clades[clade].poss.index(pos)]
                # del self._clades[clade].muts[self._clades[clade].poss.index(pos)]
                # del self._clades[clade].poss[self._clades[clade].poss.index(pos)]
        return set(self._clades[clade].poss)-set(del_set)

    def _find_closest_time(self,clade):
        # parent = self._clades[clade].parent
        if self._clades[clade].time =='':
            parent = self._clades[clade].parent
            return self._find_closest_time(parent)
        # return [clade,float(self._clades[parent].time)]
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
        # 通过突变名称从pos db中获取位点的突变类型
        # pos_list = []
        poss = []
        alts = []
        # for name in re.split(',\s*', names):
        for name in names:
            # print("/"+name+"/")
            # p0, pc = None, []
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
                    # self.posdb.update_fields(pos,['onref','name'],['*',pd.name+'&'])
                    # self.posdb.onref[pos] = '*'
                    self.posdb.add_onref(pos,'*')
                    self.posdb.add_pos_clade(pos,haplogroup)
                    # self.posdb.update_field(pos,'haplogroup',haplogroup)
                    # self.posdb[pos].haplogroup.append(haplogroup)
                else:
                    pos = nm
                    alt = '?'
                    # self.posdb.add_pos([pos,name,'?','N','*','',[haplogroup]])
                    # self.posdb.add_pos([pos,name,'?','N','*'])
                    self.posdb.add_onref(pos,'*')
                    self.posdb.add_pos([pos,name[:-1],'?','?'])
                    self.posdb.add_pos_clade(pos,haplogroup)
                    print(nm+" not found!")
                    agetext.write(nm+" not found!")
            elif name[-2] == '=':
                nm  = nms[0] if len(nms)>1 else nms[0][:-2]
                if nm in self.posdb.name_pos.keys():
                    pos = self.posdb.name_pos[nm]
                    self.posdb.add_pos_clade(pos,haplogroup)
                    self.posdb.add_onref(pos,'')
                    alt = name[-1]
                    # self.posdb[pos].haplogroup.append(haplogroup)
                else:
                    pos = nm
                    alt = name[-1]
                    # self.posdb.add_pos([pos,name,'N',alt,'','',[haplogroup]])
                    # self.posdb.add_pos([pos,name,'N',alt,''])
                    self.posdb.add_pos([pos,name[:-2],'N',alt])
                    self.posdb.add_pos_clade(pos,haplogroup)
                    self.posdb.add_onref(pos,'')
                    print(nm+" not found!")
                    agetext.write(nm+" not found!")
            else:
                nm = nms[0]
                if nm in self.posdb.name_pos.keys():
                    pos = self.posdb.name_pos[nm]
                    self.posdb.add_pos_clade(pos,haplogroup)
                    self.posdb.add_onref(pos,'')
                    pd = self.posdb[pos]
                    # p = '{}\t{}\t{}'.format(pd.pos, pd.ref, pd.alt)
                    if self.posdb.onref[pos]=='*':
                        alt = pd.ref[0]
                    else:
                        alt = pd.alt[0]
                    # self.posdb[pos].haplogroup.append(haplogroup)
                else:
                    pos =nm
                    alt = '?'
                    # self.posdb.add_pos([pos,name,'N','?','','',[haplogroup]])
                    # self.posdb.add_pos([pos,name,'N','?',''])
                    self.posdb.add_pos([pos,name,'N','?'])
                    self.posdb.add_pos_clade(pos,haplogroup)
                    self.posdb.add_onref(pos,'')
                    print(nm+" not found!")
                    agetext.write(nm+" not found!")
            # pos_list.append(p)
            poss.append(pos)
            alts.append(alt)
        return [poss,alts]
        # return set(pos_list)
            # for nm in name.split('/'):
            #     # print("/"+nm+"/")
            #     if nm[-1] == '#':
            #         break
            #     elif nm[-1] == '&':
            #         if nm[:-1] in self.posdb.name_pos.keys():
            #             pos = self.posdb.name_pos[nm[:-1]]
            #             pd = self.posdb[pos]
            #             # p= '{}\t{}\t{}'.format(pd.pos, pd.alt, pd.ref)
            #             p = '{}:{}'.format(pd.pos, pd.ref)
            #             # print(nm[:-1])
            #             # print(p)
            #             # if p0 is not None and p != p0:
            #             #     pc.append(nm+":"+p)
            #             # else:
            #             #     p0 = p
            #     elif nm[-2] == '=':
            #         if nm[:-2] in self.posdb.name_pos.keys():
            #             pos = self.posdb.name_pos[nm[:-2]]
            #             pd = self.posdb[pos]
            #             # p = '{}\t{}\t{}'.format(pd.pos, pd.ref, nm[-1])
            #             p = '{}:{}'.format(pd.pos, nm[-1])
            #             # if p0 is not None and p != p0:
            #             #     pc.append(nm+":"+p)
            #             # else:
            #             #     p0 = p
            #     else:
            #         if nm in self.posdb.name_pos.keys():
            #             pos = self.posdb.name_pos[nm]
            #             pd = self.posdb[pos]
            #             # p = '{}\t{}\t{}'.format(pd.pos, pd.ref, pd.alt)
            #             p = '{}:{}'.format(pd.pos, pd.alt)
            #             # print(nm)
            #             # print(p)
            #             # if p0 is not None and p != p0:
            #             #     pc.append(nm+":"+p)
            #             # else:
            #             #     p0 = p
        #     if p is None:
        #         pass
        #         # print('WARNING: Not find pos "{}" in posdb!'.format(name))
        #     else:
        #         # if len(pc):
        #         #     print('WARNING: Conflict pos in the mut "{}" of clade name "{}"!'.format(pc, name+":"+p0) )
        #         pos_list.append(p)
        # return set(pos_list)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="输入样本vcf，至少两个样本，所有样本与第一个样本组合计算共祖时间，需提供样本的单倍型。")
    parser.add_argument("-t", "--tree", help="输入参考tree文件路径，默认是/mnt/ilustre/users/sanger-dev/sg-users/yuguo/ytree/age/tree", required=False, default ="/mnt/ilustre/users/sanger-dev/sg-users/yuguo/ytree/age/tree" )
    parser.add_argument("-p", "--posdb", help="输入参考位点库文件路径，默认是/mnt/ilustre/users/sanger-dev/sg-users/yuguo/ytree/age/posdb", required=False, default = "/mnt/ilustre/users/sanger-dev/sg-users/yuguo/ytree/age/posdb")
    parser.add_argument("-b", "--badpos", help="输入参考坏位点库文件路径，默认是/mnt/ilustre/users/sanger-dev/sg-users/yuguo/ytree/age/badpos", required=False, default = "/mnt/ilustre/users/sanger-dev/sg-users/yuguo/ytree/age/badpos")
    parser.add_argument("-v", "--vcf", help="输入样本vcf路径，默认是当前目录下的samples-filter-8.4.vcf", required=False,default = "samples-filter-8.4.vcf")
    parser.add_argument("-l", "--bam_list", help="样本bam路径列表，--tool为build时可选，默认是当前目录下的bam.list", required=False,default = "bam.list")
    parser.add_argument("-g", "--haplogroup", help="样本所属单倍型，必须指定", required=True)
    parser.add_argument("-d", "--dep", help="样本位点过滤深度阈值，默认5", required=False, default =5)
    parser.add_argument("-o", "--outgroup", help="第几个样本开始是outgroup，默认第3个开始", required=False, default =3)
    parser.add_argument("-u", "--output", help="创建一个文件用于储存输出文件", required=True)
    args = vars(parser.parse_args())

    agetext = open(args["output"], 'w')
    print('BEGIN>>')
    agetext.write("\nBEGIN>>")
    tr = Tree(args["tree"], args['posdb'], args['badpos'],args['vcf'],args['bam_list'], args['haplogroup'], output = 'ytree',dep = args['dep'],outgroup = args['outgroup'])
    agetext.write("\n{}".format(tr))
    # tr = Tree('tree', posdb)
    print(tr)
        # posdb = PosDb('pos.db')
        # print(posdb)
        # tr = Tree(args["tree"],args['posdb'])
        # tr.locate_sample(args['sample'])

        ## test
        # tr = Tree('test/Tree.tr', 'test/pos.db')
        # print(tr)
        # # print(tr['A1'].pos)
        # # with open("sample.pos", 'w') as s:
        # #     for p in tr['O2a2b2a2a1'].pos:
        # #         s.write(p+"\n")
        #
        # print(tr.locate_sample('test/sample.pos'))
        # # print(tr.locate_sample('test/s1.pos').name)
        # tr.tree2phylip()
        # tr.tree2tnt()

    # posdb = PosDb('../posdb')
    # print(posdb)
    # with open("pos.pk",'wb') as ppfile:
    #     pickle.dump(posdb,ppfile)


    # posdb = pickle.load(open("pos.pk",'rb'))
    # vcf = Vcf('samples.vcf', posdb)
    # vcf = Vcf('samples-filter-8.4.vcf', posdb,filter_bad=True)
    # print(vcf)
    # with open("vcf.pk",'wb') as ppfile:
    #     pickle.dump(vcf,ppfile)

    # posdb = pickle.load(open("pos.pk",'rb'))
    # vcf = pickle.load(open("vcf.pk",'rb'))
    # tr = Tree('../tree', '../posdb', 'samples-filter-8.4.vcf', output = 'ytree')
    # # tr = Tree('tree', posdb)
    # print(tr)
    # tr.tree2fasta()
    # tr.run_mega()
    # with open("tree.pk",'wb') as ppfile:
    #     pickle.dump(tr,ppfile)
    # with open("pos.pk",'wb') as ppfile:
    #     pickle.dump(posdb,ppfile)

    # posdb = pickle.load(open("pos.pk",'rb'))
    # tr = pickle.load(open("tree.pk",'rb'))
    # newtr = NewTree("mega_tree.nwk", tr,calc_age=True)

    # ----------------------------------
    # tr = Tree('tree', template = True)
    # ----------------------------------
    # run_mega('tree')
    print('<<END.')
    agetext.write("\n<<END.")
    agetext.close()

        # replace_tree_name('outtree','tree.name')
