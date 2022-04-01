# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
# last_modify:20180503

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from biocluster.config import Config
import gevent
import shutil
import time
import os
import re


class GeneAnnoV2Module(Module):
    """
    用于对基因比对到nr，go，kegg，nog等注释数据库
    last modified by hongdong @ 20190215 add pfam+interpro
    """
    def __init__(self, work_id):
        super(GeneAnnoV2Module, self).__init__(work_id)
        options = [
            {"name": "fasta", "type": "string"},
            {"name": "sample_id", "type": "string"}  # 该参数可以不输入
        ]
        self.add_option(options)
        self.split_fastawgs = self.add_tool("wgs.split_fastawgs")
        self.diamond_nr = []
        self.diamond_uniport = []
        self.diamond_go = []
        self.diamond_kegg = []
        self.emapper_nog = []
        self.kegg_kobas = []
        self.pfam_scans = []
        self.interpro_scans = []
        self.nr_anno = self.add_tool("wgs.nr_anno")
        self.kegg_anno = self.add_tool("wgs.kegg_anno")
        self.go_anno = self.add_tool("wgs.go_anno")
        self.uniprot_anno = self.add_tool("wgs.uniprot_anno")
        self.eggnog_anno = self.add_tool("wgs.eggnog_anno")
        self.pfam_anno = self.add_tool("wgs.pfam_anno")
        self.interpro_anno = self.add_tool("wgs.interpro_anno")
        self.nr_db = ''
        self.uniport_db = ''
        self.go_db = ''
        self.kegg_db = ''
        self.nog_db = ''
        self.pfam_db = ''

    def check_options(self):
        if not self.option("fasta"):
            raise OptionError("缺少fasta参数")
        return True

    def split_fastawgs_run(self):
        self.split_fastawgs.set_options({
            "fasta": self.option("fasta")
        })
        self.split_fastawgs.run()

    def get_fasta_list(self):
        fasta_list = []
        with open(os.path.join(self.split_fastawgs.output_dir, "fasta.list"), 'r') as r:
            data = r.readlines()
            for line in data:
                temp = line.strip()
                if not os.path.exists(temp):
                    raise Exception({"{}文件不存在！".format(temp)})
                fasta_list.append(temp)
        return fasta_list

    def diamond_nr_run(self):
        fasta_list = self.get_fasta_list()
        for m in fasta_list:
            nr = self.add_tool("wgs.diamond")
            nr.set_options({
                "database": Config().SOFTWARE_DIR + "/" + self.nr_db.lstrip('/'),
                "query": m,
                "db_type": 'nr'
            })
            self.diamond_nr.append(nr)
        for j in range(len(self.diamond_nr)):
            self.diamond_nr[j].on("end", self.set_output, 'diamond')
        if self.diamond_nr:
            if len(self.diamond_nr) > 1:
                self.on_rely(self.diamond_nr, self.nr_anno_run)
            elif len(self.diamond_nr) == 1:
                self.diamond_nr[0].on('end', self.nr_anno_run)
        else:
            raise Exception("diamond_nr列表为空！")
        for tool in self.diamond_nr:
            gevent.sleep(1)
            tool.run()

    def diamond_kegg_run(self):
        fasta_list = self.get_fasta_list()
        for m in fasta_list:
            kegg = self.add_tool("wgs.diamond")
            kegg.set_options({
                "database": Config().SOFTWARE_DIR + "/" + self.kegg_db.lstrip('/'),
                "query": m,
                "outfmt": 6,
                "db_type": 'kegg'
            })
            self.diamond_kegg.append(kegg)
        for j in range(len(self.diamond_kegg)):
            self.diamond_kegg[j].on("end", self.set_output, 'diamond')
        if self.diamond_kegg:
            if len(self.diamond_kegg) > 1:
                self.on_rely(self.diamond_kegg, self.kegg_kobas_run)
            elif len(self.diamond_kegg) == 1:
                self.diamond_kegg[0].on('end', self.kegg_kobas_run)
        else:
            raise Exception("diamond_kegg列表为空！")
        for tool in self.diamond_kegg:
            gevent.sleep(1)
            tool.run()

    def kegg_kobas_run(self):
        """
        输入的是diamond_kegg比对结果
        :return:
        """
        self.make_blast_file('kegg')
        if len(self.diamond_kegg) > 1:
            kegg_diamond = self.output_dir + "/diamond_kegg"
        else:
            kegg_diamond = self.diamond_kegg[0].output_dir
        for m in os.listdir(kegg_diamond):
            kobas = self.add_tool("wgs.kobas_anno")
            kobas.set_options({
                "kegg_diamond_result": os.path.join(kegg_diamond, m)
            })
            self.kegg_kobas.append(kobas)
        for j in range(len(self.kegg_kobas)):
            self.kegg_kobas[j].on("end", self.set_output, 'kegg_kobas')
        if self.kegg_kobas:
            if len(self.kegg_kobas) > 1:
                self.on_rely(self.kegg_kobas, self.kegg_anno_run)
            elif len(self.kegg_kobas) == 1:
                self.kegg_kobas[0].on('end', self.kegg_anno_run)
        else:
            raise Exception("diamond_kegg列表为空！")
        for tool in self.kegg_kobas:
            gevent.sleep(1)
            tool.run()

    def diamond_go_run(self):
        fasta_list = self.get_fasta_list()
        for m in fasta_list:
            go = self.add_tool("wgs.diamond")
            go.set_options({
                "database": Config().SOFTWARE_DIR + "/" + self.go_db.lstrip('/'),
                "query": m,
                "db_type": "go"
            })
            self.diamond_go.append(go)
        for j in range(len(self.diamond_go)):
            self.diamond_go[j].on("end", self.set_output, 'diamond')
        if self.diamond_go:
            if len(self.diamond_go) > 1:
                self.on_rely(self.diamond_go, self.go_anno_run)
            elif len(self.diamond_go) == 1:
                self.diamond_go[0].on('end', self.go_anno_run)
        else:
            raise Exception("diamond_go列表为空！")
        for tool in self.diamond_go:
            gevent.sleep(1)
            tool.run()

    def diamond_uniport_run(self):
        fasta_list = self.get_fasta_list()
        for m in fasta_list:
            uniport = self.add_tool("wgs.diamond")
            uniport.set_options({
                "database": Config().SOFTWARE_DIR + '/' + self.uniport_db.lstrip('/'),
                "query": m,
                "db_type": "uniprot"
            })
            self.diamond_uniport.append(uniport)
        for j in range(len(self.diamond_uniport)):
            self.diamond_uniport[j].on("end", self.set_output, 'diamond')
        if self.diamond_uniport:
            if len(self.diamond_uniport) > 1:
                self.on_rely(self.diamond_uniport, self.uniprot_anno_run)
            elif len(self.diamond_uniport) == 1:
                self.diamond_uniport[0].on('end', self.uniprot_anno_run)
        else:
            raise Exception("diamond_uniport列表为空！")
        for tool in self.diamond_uniport:
            gevent.sleep(1)
            tool.run()

    def emapper_nog_run(self):
        fasta_list = self.get_fasta_list()
        for m in fasta_list:
            nog = self.add_tool("wgs.emapper")
            nog.set_options({
                "dmnd_db": Config().SOFTWARE_DIR + '/' + self.nog_db.lstrip('/'),
                "db_path": os.path.dirname(Config().SOFTWARE_DIR + '/' + self.nog_db.lstrip('/')),
                "fasta": m
            })
            self.emapper_nog.append(nog)
        # self.logger.info(self.emapper_nog)
        for j in range(len(self.emapper_nog)):
            self.emapper_nog[j].on("end", self.set_output, 'diamond_nog')
        if self.emapper_nog:
            if len(self.emapper_nog) > 1:
                self.on_rely(self.emapper_nog, self.eggnog_anno_run)
            elif len(self.emapper_nog) == 1:
                self.emapper_nog[0].on('end', self.eggnog_anno_run)
        else:
            raise Exception("emapper_nog列表为空！")
        for tool in self.emapper_nog:
            gevent.sleep(1)
            tool.run()

    def pfam_scan_run(self):
        fasta_list = self.get_fasta_list()
        for m in fasta_list:
            pfam_scan = self.add_tool('wgs.pfam_scan')
            pfam_scan.set_options({
                "query": m,
                'database': Config().SOFTWARE_DIR + '/' + self.pfam_db.lstrip('/')
            })
            self.pfam_scans.append(pfam_scan)
        for j in range(len(self.pfam_scans)):
            self.pfam_scans[j].on('end', self.set_output, 'pfam_scan')
        if self.pfam_scans:
            if len(self.pfam_scans) > 1:
                self.on_rely(self.pfam_scans, self.pfam_anno_run)
            elif len(self.pfam_scans) == 1:
                self.pfam_scans[0].on('end', self.pfam_anno_run)
        else:
            raise Exception("pfam_scans列表为空！")
        for tool in self.pfam_scans:
            gevent.sleep(1)
            tool.run()

    def interpro_scan_run(self):
        fasta_list = self.get_fasta_list()
        for m in fasta_list:
            pfam_scan = self.add_tool('wgs.interpro_scan')
            pfam_scan.set_options({
                "query": m
            })
            self.interpro_scans.append(pfam_scan)
        for j in range(len(self.interpro_scans)):
            self.interpro_scans[j].on('end', self.set_output, 'interpro_scan')
        if self.interpro_scans:
            if len(self.interpro_scans) > 1:
                self.on_rely(self.interpro_scans, self.interpro_anno_run)
            elif len(self.interpro_scans) == 1:
                self.interpro_scans[0].on('end', self.interpro_anno_run)
        else:
            raise Exception("interpro_scans列表为空！")
        for tool in self.interpro_scans:
            gevent.sleep(1)
            tool.run()

    def nr_anno_run(self):
        self.make_blast_file('nr')
        self.nr_anno.set_options({
            "fasta.list": os.path.join(self.split_fastawgs.output_dir, "fasta.list"),
            "nr_path": self.output_dir + "/diamond_nr"
        })
        # self.nr_anno.on('end', self.set_output, 'nr_anno')
        self.nr_anno.run()

    def kegg_anno_run(self):
        self.kegg_anno.set_options({
            "fasta.list": os.path.join(self.split_fastawgs.output_dir, "fasta.list"),
            "kegg_path": self.output_dir + "/kegg_kobas"
        })
        # self.kegg_anno.on('end', self.set_output, 'kegg_anno')
        self.kegg_anno.run()

    def go_anno_run(self):
        self.make_blast_file()
        self.go_anno.set_options({
            "fasta.list": os.path.join(self.split_fastawgs.output_dir, "fasta.list"),
            "go_path": self.output_dir + "/diamond_go"
        })
        # self.go_anno.on('end', self.set_output, 'go_anno')
        self.go_anno.run()

    def uniprot_anno_run(self):
        self.make_blast_file('uniprot')
        self.uniprot_anno.set_options({
            "fasta.list": os.path.join(self.split_fastawgs.output_dir, "fasta.list"),
            "uniprot_path": self.output_dir + "/diamond_uniport"
        })
        # self.uniprot_anno.on('end', self.set_output, 'uniprot_anno')
        self.uniprot_anno.run()

    def eggnog_anno_run(self):
        self.eggnog_anno.set_options({
            "fasta.list": os.path.join(self.split_fastawgs.output_dir, "fasta.list"),
            "eggnog_path": self.output_dir + "/diamond_nog"
        })
        # self.eggnog_anno.on('end', self.set_output, 'eggnog_anno')
        self.eggnog_anno.run()

    def pfam_anno_run(self):
        self.pfam_anno.set_options({
            "fasta.list": os.path.join(self.split_fastawgs.output_dir, "fasta.list"),
            "pfam_path": self.output_dir + "/pfam_scan"
        })
        self.pfam_anno.run()

    def interpro_anno_run(self):
        self.interpro_anno.set_options({
            "fasta.list": os.path.join(self.split_fastawgs.output_dir, "fasta.list"),
            "interpro_path": self.output_dir + "/interpro_scan"
        })
        self.interpro_anno.run()

    def linkdir(self, dirpath, dirname):
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
                # self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                # self.logger.info('cp -r %s %s' % (oldfiles[i], newdir))
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'diamond_nr':
            self.linkdir(obj.output_dir, 'diamond_nr')
        elif event['data'] == 'diamond_kegg':
            self.linkdir(obj.output_dir, 'diamond_kegg')
        elif event['data'] == "diamond_go":
            self.linkdir(obj.output_dir, 'diamond_go')
        elif event['data'] == "diamond_uniport":
            self.linkdir(obj.output_dir, 'diamond_uniport')
        elif event['data'] == "diamond_nog":
            self.linkdir(obj.output_dir, 'diamond_nog')
        elif event['data'] == "nr_anno":
            self.linkdir(obj.output_dir, 'nr_anno')
        elif event['data'] == "kegg_anno":
            self.linkdir(obj.output_dir, 'kegg_anno')
        elif event['data'] == "go_anno":
            self.linkdir(obj.output_dir, 'go_anno')
        elif event['data'] == "uniprot_anno":
            self.linkdir(obj.output_dir, 'uniprot_anno')
        elif event['data'] == "eggnog_anno":
            self.linkdir(obj.output_dir, 'eggnog_anno')
        elif event['data'] == "diamond":
            self.linkdir(obj.output_dir, 'diamond')
        elif event['data'] == "kegg_kobas":
            self.linkdir(obj.output_dir, 'kegg_kobas')
        elif event['data'] == "pfam_anno":
            self.linkdir(obj.output_dir, 'pfam_anno')
        elif event['data'] == "interpro_anno":
            self.linkdir(obj.output_dir, 'interpro_anno')
        elif event['data'] == "pfam_scan":
            self.linkdir(obj.output_dir, 'pfam_scan')
        elif event['data'] == "interpro_scan":
            self.linkdir(obj.output_dir, 'interpro_scan')
        else:
            pass

    def run(self):
        """
        :return:
        """
        super(GeneAnnoV2Module, self).run()
        self.logger.info(time.time())
        self.nr_db, self.uniport_db, self.go_db, self.kegg_db, self.nog_db, self.pfam_db =\
            self.api.api("wgs.assembly").get_anno_db()
        self.logger.info(time.time())
        self.split_fastawgs.on('end', self.diamond_nr_run)
        self.split_fastawgs.on('end', self.diamond_kegg_run)
        self.split_fastawgs.on('end', self.diamond_go_run)
        self.split_fastawgs.on('end', self.diamond_uniport_run)
        self.split_fastawgs.on('end', self.emapper_nog_run)
        self.split_fastawgs.on('end', self.pfam_scan_run)
        self.split_fastawgs.on('end', self.interpro_scan_run)
        self.on_rely([self.nr_anno, self.kegg_anno, self.go_anno, self.uniprot_anno, self.eggnog_anno,
                      self.pfam_anno, self.interpro_anno], self.end)
        self.split_fastawgs_run()

    def end(self):
        self.set_anno_file(os.path.join(self.nr_anno.output_dir, "test"),
                           os.path.join(self.uniprot_anno.output_dir, "test"),
                           os.path.join(self.kegg_anno.output_dir, "test"),
                           os.path.join(self.go_anno.output_dir, "test"),
                           os.path.join(self.eggnog_anno.output_dir, "test"),
                           os.path.join(self.pfam_anno.output_dir, "test"),
                           os.path.join(self.interpro_anno.output_dir, "test"))
        # self.rm_temp_file()
        super(GeneAnnoV2Module, self).end()

    def set_anno_file(self, nr, uniport, kegg, go, eggnog, pfam, interpro):
        file_path_ = self.output_dir + "/gene_anno"
        if not os.path.exists(file_path_):
            os.mkdir(file_path_)
        if self.option("sample_id"):
            outfile_name = os.path.join(file_path_, "{}.anno.summary".format(self.option("sample_id")))
        else:
            outfile_name = os.path.join(file_path_, "anno.summary")
        code = os.system("paste {} {} {} {} {} {} {}|cut -f 1,2,3,5,6,8,9,11,12,14,15,17,18,20,21 > {}"
                         .format(nr, uniport, kegg, go, eggnog, pfam, interpro, outfile_name))
        if code == 0:
            self.logger.info("合并注释文件成功！")
        else:
            raise Exception("合并注释文件失败！")

    def rm_temp_file(self):
        for m in os.listdir(self.output_dir):
            if m in ['diamond_go', 'diamond_nog', 'diamond_uniport', 'diamond_kegg', 'diamond_nr', 'interpro_scan',
                     'kegg_kobas', "pfam_scan"]:
                code = os.system('rm -r {}'.format(m))
                if code == 0:
                    self.logger.info("{}删除成功！")
                else:
                    self.logger.info("{}删除失败！")

    def make_blast_file(self, types='go'):
        """
        解决重运行的问题，所有的go，kegg，nr，
        :param types:
        :return:
        """
        if types == 'go':
            target_path = self.output_dir + "/diamond_go"
        elif types == 'kegg':
            target_path = self.output_dir + "/diamond_kegg"
        elif types == 'nr':
            target_path = self.output_dir + "/diamond_nr"
        else:
            target_path = self.output_dir + "/diamond_uniport"
        if not os.path.exists(target_path):
            os.mkdir(target_path)
        else:
            shutil.rmtree(target_path)
            os.mkdir(target_path)
        for m in os.listdir(self.output_dir + "/diamond"):
            if types == 'go' and re.match('.*\.go\.blast$', m):
                os.link(self.output_dir + "/diamond/{}".format(m), target_path + "/{}".format(m))
            elif types == 'kegg' and re.match('.*\.kegg\.blast$', m):
                os.link(self.output_dir + "/diamond/{}".format(m), target_path + "/{}".format(m))
            elif types == 'nr' and re.match('.*\.nr\.blast$', m):
                os.link(self.output_dir + "/diamond/{}".format(m), target_path + "/{}".format(m))
            elif types == 'uniprot' and re.match('.*\.uniprot\.blast$', m):
                os.link(self.output_dir + "/diamond/{}".format(m), target_path + "/{}".format(m))
            else:
                pass
