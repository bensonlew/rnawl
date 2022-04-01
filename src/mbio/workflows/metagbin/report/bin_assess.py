# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# __modify__ = '20190211'

from biocluster.workflow import Workflow
import os
import json
from bson import ObjectId
from mbio.packages.metagbin.common_function import link_dir


class BinAssessWorkflow(Workflow):
    """
    bin评估
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(BinAssessWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "bin_path", "type": "infile", "format": "sequence.fasta"},  # bin的文件
            {"name": 'bin_list', "type": "string"},
            {"name": 'new_bin', "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.getfasta = self.add_tool("metagbin.get_bin_fasta")
        self.amphora = self.add_tool("metagbin.amphora")
        self.amphora_anno = self.add_tool('metagbin.amphora_anno')
        self.rrna = self.add_tool("metagbin.bin_rrna")
        self.checkm = self.add_tool("metagbin.checkm_bin")
        self.step.add_steps('getfasta','checkm', 'amphora', 'rrna')
        self.remote_dir = self._sheet.output

    def run(self):
        self.rrna.on("end",self.set_db)
        self.amphora_anno.on("end",self.run_bin_rrna)
        self.amphora.on("end",self.run_amphora_summ)
        self.checkm.on("end",self.run_amphora)
        self.getfasta.on("end",self.run_checkm)
        self.run_getfasta()
        super(BinAssessWorkflow, self).run()

    def run_getfasta(self):
        opts = ({
            "bin_list": self.option("bin_list"),
            "bin_path": self.option("bin_path"),
            "new_bin": self.option("new_bin"),
        })
        self.getfasta.set_options(opts)
        self.getfasta.on('end', self.set_output, "getfasta")
        self.getfasta.run()

    def run_checkm(self):
        bin_dir = os.path.dirname(self.getfasta.option('out').prop['path'])
        opts =({
            "bin_dir": bin_dir,
        })
        self.checkm.set_options(opts)
        self.checkm.on('end', self.set_output, "checkm")
        self.checkm.run()

    def run_amphora(self):
        """
        进行物种注释
        :return:
        """
        taxon = self.get_taxon(self.checkm.output_dir + "/" + self.option('new_bin') + ".bin.summary.xls")
        binname = os.path.basename(self.getfasta.option('out').prop['path']).split('.')[0]
        opts = ({
            "bin_fa": self.getfasta.option('out'),
            "kingdom": taxon,
            "bin_name": binname,
        })
        self.amphora.set_options(opts)
        self.amphora.run()

    def run_amphora_summ(self):
        opts = {
            "amphora_dir": self.amphora.output_dir,
        }
        self.amphora_anno.set_options(opts)
        self.amphora_anno.on('end', self.set_output, "amphora")
        self.amphora_anno.run()

    def run_bin_rrna(self):
        taxon = self.get_taxon(self.checkm.output_dir + "/" + self.option('new_bin') + ".bin.summary.xls")
        self.logger.info(taxon)
        bin_fa = self.getfasta.option('out').prop['path']
        opts = ({
            "bin_fa": bin_fa,
            "kingdom":taxon,
        })
        self.rrna.set_options(opts)
        self.rrna.on('end', self.set_output, "rrna")
        self.rrna.run()

    def get_taxon(self,file):
        with open(file, 'r') as f:
            lines = f.readlines()
            bin_name = lines[1].rstrip('\n\r\t').split('\t')[-1]
            print bin_name
        return bin_name

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        bin_id = ObjectId(self.option("main_id"))
        pi_path = self.api.api("metagbin.common_api")
        mongo_key = 'bin_id,len,scf_num,lon_len,n50,aver_scf,comp,cont,strain_heter,domain'
        pi_path.add_main_detail(self.output_dir + "/" + self.option('new_bin') + ".bin.summary.xls", 'bin_detail',
                                bin_id, mongo_key, has_head=True, main_name='bi_id',main_table='bin', convert_dic={'comp':float, 'cont':float},
                                update_dic={'main_id': bin_id, 'status': 'end', 'path': self.remote_dir})
        mongo_key2 = 'bin_id,total_marker,uniq_marker,a,b,c,d,e,f'
        pi_path.add_main_detail(self.output_dir + "/" + self.option('new_bin') + ".marker.xls", 'bin_marker', bin_id,
                                mongo_key2, has_head=True, main_name='bi_id')
        mongo_key3 = 'bin_id,taxon'
        pi_path.add_main_detail(self.output_dir + '/summary.anno.xls', 'bin_taxon', bin_id,
                                mongo_key3, has_head=True, main_name='bi_id')
        if os.path.exists(self.output_dir + '/all_bins.16s.xls'):
            mongo_key4 = 'sca_id,bin_id,start,end,strand,seq'
            pi_path.add_main_detail(self.output_dir + '/all_bins.16s.xls', 'bin_s16', bin_id,
                                    mongo_key4, has_head=True, main_name='bi_id')
        self.end()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'getfasta':
            for i in [self.option('new_bin') + '.fa']:
                if os.path.exists(self.output_dir + '/' + i):
                    os.remove(self.output_dir + '/' + i)
                os.link(obj.output_dir + '/' + i, self.output_dir + '/' + i)
        if event['data'] == 'checkm':
            for i in [self.option('new_bin') + ".marker.xls", self.option('new_bin') + ".bin.summary.xls"]:
                if os.path.exists(self.output_dir + '/' + i):
                    os.remove(self.output_dir + '/' + i)
                os.link(obj.output_dir + '/' + i, self.output_dir + '/' + i)
        if event['data'] == 'amphora':
            for i in ['summary.anno.xls']:
                if os.path.exists(self.output_dir + '/' + i):
                    os.remove(self.output_dir + '/' + i)
                os.link(obj.output_dir + '/' + i, self.output_dir + '/' + i)
        if event['data'] == 'rrna':
            if os.path.exists(obj.output_dir + '/all_bins.16s.xls'):
                for i in [self.option('new_bin') + '.rRNA.gff', self.option('new_bin') + '.rRNA.fnn', 'all_bins.16s.xls']:
                    if os.path.exists(self.output_dir + '/' + i):
                        os.remove(self.output_dir + '/' + i)
                    os.link(obj.output_dir + '/' + i, self.output_dir + '/' + i)

    def end(self):
        bin_name = self.option("new_bin")
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "校正后%s结果目录" % bin_name],
            ["*.bin.summary.xls", "", "%s的基本信息统计表" % bin_name],
            ["*.fa", "", "校正后bin的序列文件"],
            ["summary.anno.xls", "", "校正后%s注释物种分类文件" % bin_name],
            ["*.marker.xls", "", "校正后%s的marker评估表" % bin_name],
        ])
        super(BinAssessWorkflow, self).end()