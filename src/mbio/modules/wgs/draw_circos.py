# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
# last_modify:20180414

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from biocluster.config import Config
from bson.objectid import ObjectId
import gevent
import re
import os


class DrawCircosModule(Module):
    """
    用于对snp indel vcf文件进行注释统计
    """
    def __init__(self, work_id):
        super(DrawCircosModule, self).__init__(work_id)
        options = [
            {"name": "windows", "type": "int", "default": 100000},
            {"name": "snp", "type": "string"},  # snp.anno.primary.vcf
            {"name": "indel", "type": "string"},  # indel.anno.primary.vcf
            {"name": "genome_version_id", "type": "string"},
            {"name": "sv_path", "type": "string"},
            {"name": "cnv_path", "type": "string"},
            # {"name": "chrlist", "type": "string"},
            {"name": "chrs", "type": "string", "default": "all"}  # 传进来是染色体列，如：chr1,chr2,chr3 逗号分隔
        ]
        self.add_option(options)
        self.cnv_sv = {}
        self.circos_tools = []
        self.total_chrlist = ''
        self.snpeff_path = ''

    def check_options(self):
        if not self.option('snp'):
            raise OptionError('必须提供snp结果表', code="24500601")
        if not self.option('indel'):
            raise OptionError('必须提供indel结果表', code="24500602")
        # if not self.option("gff"):
        #     raise OptionError("必须提供gff结果文件")
        # if not self.option("chrlist"):
        #     raise OptionError("必须提供chrlist文件！")
        if not self.option("sv_path"):
            raise OptionError("必须提供sv_path文件！", code="24500603")
        if not self.option("cnv_path"):
            raise OptionError("必须提供cnv_path文件！", code="24500604")
        # if not self.option("chrs"):
        #     raise OptionError("必须提供chrs文件！")
        return True

    def new_chrlist(self):
        self.logger.info("开始设置chrlist文件！")
        if self.option("chrs").lower() == 'all':
            chrs = self.get_chrs()
        else:
            chrs = self.option("chrs")
        n = 0
        with open(self.total_chrlist, 'r') as r, open(self.work_dir + '/chrlist.txt', "w") as w:
            data = r.readlines()
            for line in data:
                # self.logger.info(line)
                tmp = line.strip().split('\t')
                if tmp[0].strip() in chrs.strip().split(','):
                    n += 1
                    # self.logger.info(n)
                    w.write(line)
        if n != len(chrs.split(',')) or n == 0:
            self.set_error("页面上面选择的染色体在chrlist中不存在！", code="24500611")
        self.logger.info("设置chrlist文件成功！")

    def get_chrs(self):
        """
        如果前端传进来是all 默认计算所有的染色体，没有染色体的时候，计算前二十个scaflods
        :return:
        """
        chrs = []
        scas = []
        with open(self.total_chrlist, 'r') as r:
            data = r.readlines()
            for line in data:
                temp = line.strip().split("\t")
                if re.match(r'^chr.*', temp[0]):
                    chrs.append(temp[0])
                else:
                    scas.append(temp[0])
        if len(chrs) != 0:
            result = ','.join(chrs)
        else:
            result = ','.join(scas[0:20])
        self.logger.info(result[0:20])
        return result

    def set_cnv_sv_file(self):
        """
        {"A8_10": {"cnv": "", "sv": ""}}
        :return:
        """
        self.logger.info("开始设置cnv_sv文件！")
        for file_ in os.listdir(self.option("cnv_path")):
            m = re.match(r'(.*)\.cnv\.anno\.xls$', file_)
            if m:
                self.cnv_sv[m.group(1)] = {}
                self.cnv_sv[m.group(1)]['cnv'] = os.path.join(self.option("cnv_path"), file_)
                self.api.api("wgs.api_base").check_exists(os.path.join(self.option("sv_path"),
                                                                       "{}.sv.anno.xls".format(m.group(1))))
                self.cnv_sv[m.group(1)]['sv'] = os.path.join(self.option("sv_path"), "{}.sv.anno.xls".format(m.group(1)))
        self.logger.info(self.cnv_sv)
        self.logger.info("设置cnv_sv文件成功！")

    def draw_circos(self):
        self.new_chrlist()
        self.set_cnv_sv_file()
        for key in self.cnv_sv.keys():
            circos = self.add_tool("wgs.draw_circos")
            circos.set_options({
                "windows": self.option("windows"),
                "snp": self.option("snp"),
                "indel": self.option("indel"),
                "gff": self.snpeff_path,
                "sv": self.cnv_sv[key]['sv'],
                "chrlist": self.work_dir + '/chrlist.txt',
                "cnv": self.cnv_sv[key]['cnv']
            })
            self.circos_tools.append(circos)
        for j in range(len(self.circos_tools)):
            self.circos_tools[j].on("end", self.set_output, 'circos')
        if self.circos_tools:
            if len(self.circos_tools) > 1:
                self.on_rely(self.circos_tools, self.end)
            elif len(self.circos_tools) == 1:
                self.circos_tools[0].on('end', self.end)
        else:
            self.set_error("circos_tools列表为空！", code="24500612")
        for tool in self.circos_tools:
            gevent.sleep(1)
            tool.run()

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
        if event['data'] == 'circos':
            self.linkdir(obj.output_dir, self.output_dir)
        else:
            pass

    def run(self):
        super(DrawCircosModule, self).run()
        self.logger.info("开始这是gff与chrlist文件")
        self.get_gene_eff()
        self.logger.info("设置gff与chrlist文件成功")
        self.draw_circos()

    def get_gene_eff(self):
        result = self.api.api('wgs.api_base').col_find_one("sg_species_version",
                                                           {"_id": ObjectId(self.option("genome_version_id"))})
        if result:
            self.snpeff_path = Config().SOFTWARE_DIR + '/database/dna_wgs_geneome/' + \
                               os.path.join(os.path.dirname(result['snpeff_path']), "ref.gff")
            if not os.path.exists(self.snpeff_path):
                self.set_error("文件%s不存在!", variables=(self.snpeff_path), code="24500613")
            self.total_chrlist = Config().SOFTWARE_DIR + '/database/dna_wgs_geneome/' + result['total_chrlist']
            if not os.path.exists(self.total_chrlist):
                self.set_error("文件%s不存在!", variables=(self.total_chrlist), code="24500614")
        else:
            self.set_error("在sg_species_version中没有找到%s对应信息", variables=(self.option("genome_version_id")), code="24500615")

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(DrawCircosModule, self).end()
