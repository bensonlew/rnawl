# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
# last_modify:20180503

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from biocluster.config import Config
import gevent
import time
import os


class GeneAnnoModule(Module):
    """
    用于对基因比对到nr，go，kegg，nog等注释数据库
    """
    def __init__(self, work_id):
        super(GeneAnnoModule, self).__init__(work_id)
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
        self.nr_anno = self.add_tool("wgs.nr_anno")
        self.kegg_anno = self.add_tool("wgs.kegg_anno")
        self.go_anno = self.add_tool("wgs.go_anno")
        self.uniprot_anno = self.add_tool("wgs.uniprot_anno")
        self.eggnog_anno = self.add_tool("wgs.eggnog_anno")
        self.nr_db = ''
        self.uniport_db = ''
        self.go_db = ''
        self.kegg_db = ''
        self.nog_db = ''

    def check_options(self):
        if not self.option("fasta"):
            raise OptionError("缺少fasta参数", code="24500701")
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
                    self.set_error({"%s文件不存在！", variables=(temp)}, code="24500715")
                fasta_list.append(temp)
        return fasta_list

    def diamond_nr_run(self):
        fasta_list = self.get_fasta_list()
        for m in fasta_list:
            nr = self.add_tool("wgs.diamond")
            nr.set_options({
                "database": Config().SOFTWARE_DIR + "/" + self.nr_db.lstrip('/'),
                "query": m
            })
            self.diamond_nr.append(nr)
        for j in range(len(self.diamond_nr)):
            self.diamond_nr[j].on("end", self.set_output, 'diamond_nr')
        if self.diamond_nr:
            if len(self.diamond_nr) > 1:
                self.on_rely(self.diamond_nr, self.nr_anno_run)
            elif len(self.diamond_nr) == 1:
                self.diamond_nr[0].on('end', self.nr_anno_run)
        else:
            self.set_error("diamond_nr列表为空！", code="24500716")
        for tool in self.diamond_nr:
            gevent.sleep(1)
            tool.run()

    def diamond_kegg_run(self):
        fasta_list = self.get_fasta_list()
        for m in fasta_list:
            kegg = self.add_tool("wgs.diamond")
            kegg.set_options({
                "database": Config().SOFTWARE_DIR + "/" + self.kegg_db.lstrip('/'),
                "query": m
            })
            self.diamond_kegg.append(kegg)
        for j in range(len(self.diamond_kegg)):
            self.diamond_kegg[j].on("end", self.set_output, 'diamond_kegg')
        if self.diamond_kegg:
            if len(self.diamond_kegg) > 1:
                self.on_rely(self.diamond_kegg, self.kegg_anno_run)
            elif len(self.diamond_kegg) == 1:
                self.diamond_kegg[0].on('end', self.kegg_anno_run)
        else:
            self.set_error("diamond_kegg列表为空！", code="24500717")
        for tool in self.diamond_kegg:
            gevent.sleep(1)
            tool.run()

    def diamond_go_run(self):
        fasta_list = self.get_fasta_list()
        for m in fasta_list:
            go = self.add_tool("wgs.diamond")
            go.set_options({
                "database": Config().SOFTWARE_DIR + "/" + self.go_db.lstrip('/'),
                "query": m
            })
            self.diamond_go.append(go)
        for j in range(len(self.diamond_go)):
            self.diamond_go[j].on("end", self.set_output, 'diamond_go')
        if self.diamond_go:
            if len(self.diamond_go) > 1:
                self.on_rely(self.diamond_go, self.go_anno_run)
            elif len(self.diamond_go) == 1:
                self.diamond_go[0].on('end', self.go_anno_run)
        else:
            self.set_error("diamond_go列表为空！", code="24500718")
        for tool in self.diamond_go:
            gevent.sleep(1)
            tool.run()

    def diamond_uniport_run(self):
        fasta_list = self.get_fasta_list()
        for m in fasta_list:
            uniport = self.add_tool("wgs.diamond")
            uniport.set_options({
                "database": Config().SOFTWARE_DIR + '/' + self.uniport_db.lstrip('/'),
                "query": m
            })
            self.diamond_uniport.append(uniport)
        for j in range(len(self.diamond_uniport)):
            self.diamond_uniport[j].on("end", self.set_output, 'diamond_uniport')
        if self.diamond_uniport:
            if len(self.diamond_uniport) > 1:
                self.on_rely(self.diamond_uniport, self.uniprot_anno_run)
            elif len(self.diamond_uniport) == 1:
                self.diamond_uniport[0].on('end', self.uniprot_anno_run)
        else:
            self.set_error("diamond_uniport列表为空！", code="24500719")
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
        self.logger.info(self.emapper_nog)
        for j in range(len(self.emapper_nog)):
            self.emapper_nog[j].on("end", self.set_output, 'diamond_nog')
        if self.emapper_nog:
            if len(self.emapper_nog) > 1:
                self.on_rely(self.emapper_nog, self.eggnog_anno_run)
            elif len(self.emapper_nog) == 1:
                self.emapper_nog[0].on('end', self.eggnog_anno_run)
        else:
            self.set_error("emapper_nog列表为空！", code="24500720")
        for tool in self.emapper_nog:
            gevent.sleep(1)
            tool.run()

    def nr_anno_run(self):
        self.nr_anno.set_options({
            "fasta.list": os.path.join(self.split_fastawgs.output_dir, "fasta.list"),
            "nr_path": self.output_dir + "/diamond_nr"
        })
        # self.nr_anno.on('end', self.set_output, 'nr_anno')
        self.nr_anno.run()

    def kegg_anno_run(self):
        self.kegg_anno.set_options({
            "fasta.list": os.path.join(self.split_fastawgs.output_dir, "fasta.list"),
            "kegg_path": self.output_dir + "/diamond_kegg"
        })
        # self.kegg_anno.on('end', self.set_output, 'kegg_anno')
        self.kegg_anno.run()

    def go_anno_run(self):
        self.go_anno.set_options({
            "fasta.list": os.path.join(self.split_fastawgs.output_dir, "fasta.list"),
            "go_path": self.output_dir + "/diamond_go"
        })
        # self.go_anno.on('end', self.set_output, 'go_anno')
        self.go_anno.run()

    def uniprot_anno_run(self):
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
        else:
            pass

    def run(self):
        """
        :return:
        """
        super(GeneAnnoModule, self).run()
        self.logger.info(time.time())
        self.nr_db, self.uniport_db, self.go_db, self.kegg_db, self.nog_db = self.api.api("wgs.assembly").get_anno_db()
        self.logger.info(time.time())
        self.split_fastawgs.on('end', self.diamond_nr_run)
        self.split_fastawgs.on('end', self.diamond_kegg_run)
        self.split_fastawgs.on('end', self.diamond_go_run)
        self.split_fastawgs.on('end', self.diamond_uniport_run)
        self.split_fastawgs.on('end', self.emapper_nog_run)
        self.on_rely([self.nr_anno, self.kegg_anno, self.go_anno, self.uniprot_anno, self.eggnog_anno], self.end)
        self.split_fastawgs_run()

    def end(self):
        self.set_anno_file(os.path.join(self.nr_anno.output_dir, "test"),
                           os.path.join(self.uniprot_anno.output_dir, "test"),
                           os.path.join(self.kegg_anno.output_dir, "test"),
                           os.path.join(self.go_anno.output_dir, "test"),
                           os.path.join(self.eggnog_anno.output_dir, "test"))
        self.rm_temp_file()
        super(GeneAnnoModule, self).end()

    def set_anno_file(self, nr, uniport, kegg, go, eggnog):
        file_path_ = self.output_dir + "/gene_anno"
        if not os.path.exists(file_path_):
            os.mkdir(file_path_)
        if self.option("sample_id"):
            outfile_name = os.path.join(file_path_, "{}.anno.summary".format(self.option("sample_id")))
        else:
            outfile_name = os.path.join(file_path_, "anno.summary")
        code = os.system("paste {} {} {} {} {}|cut -f 1,2,3,5,6,8,9,11,12,14,15 > {}"
                         .format(nr, uniport, kegg, go, eggnog, outfile_name))
        if code == 0:
            self.logger.info("合并注释文件成功！")
        else:
            self.set_error("合并注释文件失败！", code="24500721")

    def rm_temp_file(self):
        for m in os.listdir(self.output_dir):
            if m in ['diamond_go', 'diamond_nog', 'diamond_uniport', 'diamond_kegg', 'diamond_nr']:
                code = os.system('rm -r {}'.format(m))
                if code == 0:
                    self.logger.info("{}删除成功！")
                else:
                    self.logger.info("{}删除失败！")
