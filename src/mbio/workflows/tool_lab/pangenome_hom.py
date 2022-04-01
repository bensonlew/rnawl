# -*- coding: utf-8 -*-


from biocluster.workflow import Workflow
import os,re
import shutil
from biocluster.core.exceptions import OptionError
from Bio import SeqIO

class PangenomeHomWorkflow(Workflow):
    """

    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PangenomeHomWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "genome_dir", "type": "infile", "format": "tool_lab.input_dir"},  # 输入faa文件夹下包含蛋白序列
            {"name": "identity", "type": "float", "default": 0.005},  ##给出cdhit的参数identity
            {"name": "coverage", "type": "float", "default": 0.00001},  # 给出cdhit的参数coverage
            {"name": "method", "type": "string", "default":"OMCL"}, # 三种方法BDBH、OMCL、OCOG
            {"name": "inflation", "type": "float", "default": 1.5},  # Markov Inflation Index
            {"name": "test", "type": "string", "default": "true"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())

    def check_file(self,file, type):
        """
        对序列文件做检查
        :param file:
        :return:
        """
        for uniprot_iterator in SeqIO.parse(file, "fasta"):
            dd = os.path.basename(file)
            if type in ['fna']:
                if re.search("[^atgcnATGCN]",str(uniprot_iterator.seq)):
                    raise OptionError("核苷酸文件{}的{}含有特殊字符!".format(dd, uniprot_iterator.id))
            elif type in ['faa']:
                if re.search("[^GAVLIFWYDHNEKQMRSTCPUOX]",str(uniprot_iterator.seq)):
                    raise OptionError("氨基酸文件{}的{}含有特殊字符!".format(dd, uniprot_iterator.id))

    def run(self):
        group_table = self.work_dir+'/group.xls'
        with open(self.option('genome_dir').path+'/list.txt') as fr,  open(group_table,'w') as fw:
            fw.write('#sample\tgroup\n')
            data = fr.readlines()
            all_file = {}
            for line in data[1:]:
                sp = line.strip().split('\t')
                fw.write(sp[0]+'\t'+sp[0]+'\n')
                all_file[sp[0]] =sp[1]
        if os.path.exists(self.work_dir + "/faa"):
            shutil.rmtree(self.work_dir + "/faa")
        os.mkdir(self.work_dir + "/faa")
        for file in all_file.keys():
            if os.path.exists(self.option("genome_dir").path + "/" +all_file[file]):
                if all_file[file].endswith(".faa"):
                    self.check_file(self.option("genome_dir").path + "/" +all_file[file], "faa")
                    os.link(self.option("genome_dir").path + "/" +all_file[file],self.work_dir + "/faa/" + file+".faa")
                elif all_file[file].endswith(".fna"):
                    self.check_file(self.option("genome_dir").path + "/" + all_file[file], "fna")
                    os.link(self.option("genome_dir").path + "/" + all_file[file],
                            self.work_dir + "/faa/" + file + ".fna")
            else:
                raise OptionError("文件%s不存在".format(file))
        self.anno = self.add_module("tool_lab.pan")
        options = {
            "fasta_dir": self.work_dir + "/faa",
            "method": "get_homologus",
            "homologus_mehod": self.option("method"),
            "identity": self.option("identity"),
            "coverage": self.option("coverage"),
            "inflation": self.option("inflation"),
            "group_table": group_table
        }
        self.anno.set_options(options)
        self.anno.on('end',self.set_db)
        self.anno.run()
        super(PangenomeHomWorkflow, self).run()

    def set_db(self):
        """
        导表程序
        """
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
        shutil.copytree(self.anno.output_dir, self.output_dir)

        pangenomes_dir = os.path.join(self.output_dir, "pangenomes")
        cluster_dir = os.path.join(pangenomes_dir, "cluster")
        api_pan = self.api.api('tool_lab.pan')

        api_venn = self.api.api('tool_lab.pan_venn')
        cluster_data = os.path.join(cluster_dir, 'homologues_cluster.xls')

        self.logger.info("正在导泛基因组分析的主表")
        pan_id = self.option("main_id")
        self.logger.info("正在导泛基因组分析的cluster结果表")
        api_pan.add_pan_detail(pan_id, cluster_data,collection_type='hom')

        pan_graph = os.path.join(cluster_dir, 'pan_cluster_genome.xls')
        pan_detail = os.path.join(cluster_dir, 'pan_cluster.xls')
        core_graph = os.path.join(cluster_dir, 'core_cluster_genome.xls')
        core_detail = os.path.join(cluster_dir, 'core_cluster.xls')
        new_graph = os.path.join(cluster_dir, 'new_gene_cluster_genome.xls')
        new_detail = os.path.join(cluster_dir, 'new_gene_cluster.xls')
        self.logger.info("正在导公式的结果表")
        #先导pan再导core
        api_pan.add_pan_genome(pan_id, 'pan', pan_graph, pan_detail,collection_type='hom')
        api_pan.add_pan_genome(pan_id, 'core', core_graph, core_detail,collection_type='hom')
        api_pan.add_pan_genome(pan_id, 'newgene', new_graph, new_detail,collection_type='hom')


        self.logger.info("正在进行venn图分析")
        venn_path = self.anno.output_dir + '/pangenomes/venn/new_venn_graph.xls'
        api_venn.add_venn_detail_tool(pan_id, venn_path,collection_type='hom')
        remote_png = self.sheet.output
        remote_png_path =  {}

        if not remote_png:
            remote_png = './'
        if os.path.exists(self.output_dir+'/pangenomes/cluster/core_cluster.xls_core_Tettelin.png'):
            core_png = remote_png + '/pangenomes/cluster/core_cluster.xls_core_Tettelin.png'
            remote_png_path['core'] = core_png
        if os.path.exists(self.output_dir+'/pangenomes/cluster/pan_cluster.xls_pan.png'):
            pan_png = remote_png +'/pangenomes/cluster/pan_cluster.xls_pan.png'
            remote_png_path['pan'] = pan_png
        if os.path.exists(self.output_dir + '/pangenomes/cluster/new_gene_cluster.xls_newgene.png'):
            new_png = remote_png + '/pangenomes/cluster/new_gene_cluster.xls_newgene.png'
            remote_png_path['new'] = new_png
        api_pan.update_main(pan_id,cluster_data, collection_type='hom', png_info=remote_png_path)
        self.end()

    def end(self):


        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "dir", "所有样品基因组的泛基因组分析目录", 0, ""],
            ["./venn", "dir", "所有样品基因组的泛基因组分析venn目录", 0],
            ["./venn/venn_table.xls", "xls", "泛基因组分析venn结果", 0],
            ["./cluster", "dir", "所有样品基因组的泛基因组分析聚类同源蛋白目录", 0],
            ["./cluster/pan_cluster_genome.xls", "xls", "泛基因组分析的GET_HOMOLOGOUS的pangenomes计算公式结果表", 0],
            ["./cluster/pan_cluster.xls", "xls", "泛基因组分析的pangenomes曲线结果表", 0],
            ["./cluster/pan_cluster.xls_pan.png", "png", "泛基因组分析的pangenomes曲线结果图", 0],
            ["./cluster/new_gene_cluster_genome.xls", "xls", "newgene计算公式结果表", 0],
            ["./cluster/new_gene_cluster.xls", "xls", "newgene曲线结果表", 0],
            ["./cluster/new_gene_cluster.xls_newgene.png", "png", "newgene曲线结果图", 0],
            ["./cluster/homologues_cluster.xls", "xls", "同源基因结果表", 0],
            ["./cluster/core_cluster_genome.xls", "xls", "泛基因组分析的GET_HOMOLOGOUS的coregene计算公式结果表", 0],
            ["./cluster/core_cluster.xls", "xls", "泛基因组分析的coregene曲线结果表", 0],
            ["./cluster/core_cluster.xls_core_Tettelin.png", "png", "泛基因组分析的coregene曲线结果图", 0],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(PangenomeHomWorkflow, self).end()

if __name__ == "__main__":
    from biocluster.wsheet import Sheet
    data = {
        'name': 'tool_lab.pangenome_hom',
        'id': 'tsg_36964',
        'type': 'workflow',
        'options': {
            "main_id" : "5e4ce5e817b2bf4b326f4b20",
            "genome_dir" : "/mnt/ilustre/users/sanger-dev/sg-users/zouguanqing/small_tool/gene_fa",
            "test" : "true"
        }
    }

    wsheet = Sheet(data=data)

    wf = PangenomeHomWorkflow(wsheet)
    wf.run()
