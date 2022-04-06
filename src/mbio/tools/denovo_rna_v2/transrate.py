# coding=utf-8
import os
import glob
import random
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
import shutil

__author__ = 'liubinxu'


class TransrateAgent(Agent):
    """
    评估、过滤转录组组装结果
    """
    def __init__(self, parent):
        super(TransrateAgent, self).__init__(parent)
        options = [
            {'type': 'infile', 'name': 'assembly', 'format': 'denovo_rna_v2.trinity_fasta'},
            {'type': 'infile', 'name': 'left', 'format': 'sequence.fastq'},
            {'type': 'infile', 'name': 'right', 'format': 'sequence.fastq'},
            {'default': 20, 'type': 'int', 'name': 'threads'},
            {'type': 'bool', 'name': 'filter', 'default': False},
            {'type': 'outfile', 'name': 'exp_result', "format": "denovo_rna_v2.common"},
            {'type': 'outfile', 'name': 'good_fa', 'format': 'denovo_rna_v2.trinity_fasta'},
            {'type': 'outfile', 'name': 'bad_fa', "format": "denovo_rna_v2.common"},
            {'type': 'outfile', 'name': 'result', "format": "denovo_rna_v2.common"},
        ]
        self.add_option(options)
        self._max_rerun_time = 3
        self._memory_increase_step = 60

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = self.option('threads')
        self._memory = "{}G".format('40')

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ["assemblies.csv", "csv", "Transrate 评估结果"],
            ["good.Trinity.fasta", "fasta", "Transrate 过滤后得分较高的序列文件"],
            ["bad.Trinity.fasta", "fasta", "Transrate 过滤后得分较低的序列文件"],
            ["Trinity.fasta_quant.sf", "txt", "Transrate 统计的基因表达量文件"]
            ])
        super(TransrateAgent, self).end()



class TransrateTool(Tool):
    """
    评估、过滤转录组组装结果
    """
    def __init__(self, config):
        super(TransrateTool, self).__init__(config)
        self.software_dir = self.config.SOFTWARE_DIR
        self.hisat_path = 'bioinfo/align/hisat2/hisat2-2.1.0/'
        self.samtools_path =  '/bioinfo/align/samtools-1.3.1/'
        self.python_path = 'miniconda2/bin/python'
        self.cmd_run = '/bioinfo/denovo_rna_v2/run_cmd.sh'
        self.transrate = '/bioinfo/denovo_rna_v2/transrate-1.0.3-linux-x86_64_new/transrate'
        python_path = self.config.SOFTWARE_DIR + '/program/Python/bin/'
        self.set_environ(PATH=python_path)


    def fastq_samples(self):
        """
        抽样评估,目前数据大于30G在进行抽样
        """
        left_fq_size = os.path.getsize(self.option('left').prop['path'])
        right_fq_size = os.path.getsize(self.option('right').prop['path'])
        random_range = 1
        if left_fq_size + right_fq_size > 30 * 1024 * 1024 * 1024:
            random_range = int((left_fq_size + right_fq_size)/(30 * 1024 * 1024 *1024)) + 1
        with open(self.option('left').prop['path'], 'r') as left, \
             open(self.option('right').prop['path'], 'r') as right, \
             open('filter.right.fq', 'w') as left_w, \
             open('filter.left.fq', 'w') as right_w:
            while True:
                l1 = left.readline()
                if l1:
                    pass
                else:
                    break
                l2 = left.readline()
                l3 = left.readline()
                l4 = left.readline()
                r1 = right.readline()
                r2 = right.readline()
                r3 = right.readline()
                r4 = right.readline()
                r = random.randint(1, random_range)
                if r == 1:
                    left_w.write(l1)
                    left_w.write(l2)
                    left_w.write(l3)
                    left_w.write(l4)
                    right_w.write(r1)
                    right_w.write(r2)
                    right_w.write(r3)
                    right_w.write(r4)
                else:
                    pass

    def hisat_build(self):
        """
        建立索引
        """
        cmd = "{}hisat2-build -p {} -f {} ref_index".format(
            self.hisat_path,
            self.option('threads'),
            self.option("assembly").prop['path'])
        self.logger.info("开始运行hisat2-build，进行建索引")
        command = self.add_command("hisat_build", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("hisat-build运行完成")
            os.system("touch  {}.finished".format("hisat_build"))
            return True
        else:
            command.rerun()
            self.wait(command)
            if command.return_code == 0:
                return True
            else:
                raise Exception("建立索引出错", code = "32005901")
        os.system("touch  {}.finished".format("hisat_build"))

    def hisat_mapping(self):
        """
        比对
        """
        cmd = ""
        if self.option("right").is_set:
            cmd = "{}hisat2 --no-spliced-alignment --no-mixed -p {} -q -x {} -1 {} -2 {} -S hisat_mapping.unsorted.sam".format(
                self.hisat_path,
                self.option('threads'),
                "ref_index",
                self.option("left").prop["path"],
                self.option("right").prop["path"])
        else:
            cmd = "{}hisat2 --no-spliced-alignment -p {} -q -x {} -U {} -S hisat_mapping.unsorted.sam".format(
                self.hisat_path,
                self.option('threads'),
                "ref_index",
                self.option("left").prop["path"])

        command = self.add_command("hisat_mapping", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("hisat运行完成")
            os.system("touch  {}.finished".format("hisat_mapping"))
            return True
        elif command.return_code == None:
            command.rerun()
            self.wait()
            if command.return_code == 0:
                self.logger.info("hisat运行完成")
        else:
            self.logger.error("hisat运行出错")
            self.set_error("hisat运行出错", code = "32005902")
            self.set_error("运行hisat出错", code="32005904")
        os.system("touch  {}.finished".format("hisat_mapping"))


    def sam_bam(self):
        """
        运行samtools，将hisat比对的结果文件sam文件转成bam文件
        """
        sam_path = os.path.join(self.work_dir, "hisat_mapping.unsorted.sam")
        cmd = "{}samtools view -@ {} -bS {} -o hisat_mapping.unsorted.bam".format(self.samtools_path, self.option('threads'), sam_path)
        command = self.add_command("sam2bam", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("sam2bam运行完成")
            os.system("touch  {}.finished".format("sam2bam"))
            return True
        elif command.return_code == None:
            command.rerun()
            if command.return_code == 0:
                self.logger.info("sam2bam运行完成")
        else:
            self.logger.error("sam2bam运行出错")
            self.set_error("sam2bam运行出错", code = "32005903")
        os.system("touch  {}.finished".format("sam2bam"))

    def run_transrate(self):
        bam_path = os.path.join(self.work_dir, "hisat_mapping.unsorted.bam")
        right_read = ""
        if self.option('right').is_set:
            right_read = self.option('right').prop['path']
        else:
            os.system("touch right.fq")
            right_read = "right.fq"

        cmd = '{} --assembly={} --left={} --right={} --threads={} --bam={}'.format(
            self.transrate,
            self.option('assembly').prop['path'],
            self.option('left').prop['path'],
            right_read,
            self.option('threads'),
            bam_path
        )
        if self.option('filter'):
            self.fastq_samples()
            cmd = '{} --assembly={} --left={} --right={} --threads={} --bam={}'.format(
                self.transrate,
                self.option('assembly').prop['path'],
                "filter.left.fq",
                "filter.right.fq",
                self.option('threads'),
                bam_path
            )
        else:
            pass
        cmd_name = 'transrate'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} Finished successfully".format(cmd_name))
        elif command.return_code == 1:  # add memory limit by shicaiping at 20180816
            self.add_state("memory_limit", "memory is low!")
        else:
            command.rerun()
            if command.return_code == 0:
                self.logger.info("transrate运行完成")
        os.system("touch  {}.finished".format("transrate"))

        '''
        # 不跳过slurm限制会导致机器无响应
        if command.return_code == 0:
            self.logger.info("{} Finished successfully".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("{} Failed and returned None, we will try it again.".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} Finished successfully".format(cmd_name))
            else:
                self.set_error("{} Failed. {}".format(cmd_name, cmd))
        else:
            self.logger.info("Failed try skip slurm method")
            if os.path.exists('transrate_results'):
                shutil.rmtree('transrate_results')
            else:
                pass
            with open('skip_slurm.sh', 'w') as skip_slurm:
                skip_slurm.write("ssh localhost '{} --assembly={} --left={} --right={} --threads={} --output={}'".format(
                    self.software_dir + self.transrate,
                    self.option('assembly').prop['path'],
                    self.option ('left').prop['path'],
                    self.option('right').prop['path'],
                    self.option('threads'),
                    self.work_dir + '/transrate_results'
                ))
            cmd_skip_slurm = '{} skip_slurm.sh'.format(self.cmd_run)
            command = self.add_command('cmd_skip_slurm', cmd_skip_slurm)
            command.run()
            self.wait()
            if command.return_code == 0:
                self.logger.info("{} Finished successfully".format(cmd_name))
            else:
                self.set_error("{} Failed. {}".format(cmd_name, cmd))
        '''

    def set_output(self):
        base_fa = os.path.basename(self.option('assembly').prop['path'])
        base = base_fa
        if base.endswith('.fa'):
            base = base[:-3]
        elif base.endswith('.fasta'):
            base = base[:-6]
        else:
            base = os.path.splitext(base)[0]
            # pass
        stat_files = glob.glob('transrate_results/assemblies.csv')
        good_fa = glob.glob('transrate_results/' + base + '/good.' + base_fa)
        bad_fa = glob.glob('transrate_results/' + base + '/bad.' + base_fa)
        quant_sf = glob.glob('transrate_results/' + base + '/' + base_fa + '_quant.sf')
        all_files = stat_files + good_fa + bad_fa + quant_sf
        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)
        if os.path.exists('transrate_results/' + base + '/'+ 'quant.sf'):
            link = os.path.join(self.output_dir, base_fa + '_quant.sf')
            if os.path.exists(link):
                os.remove(link)
            os.link('transrate_results/' + base + '/'+ 'quant.sf', link)
        self.option('exp_result', os.path.join(self.output_dir, base_fa + '_quant.sf'))
        self.option('good_fa', os.path.join(self.output_dir, 'good.' + base_fa))
        self.option('bad_fa', os.path.join(self.output_dir, 'bad.' + base_fa))
        self.option('result', os.path.join(self.output_dir, 'assemblies.csv'))

    def run(self):
        super(TransrateTool, self).run()
        if os.path.exists("hisat_build.finished"):
            pass
        else:
            self.hisat_build()
        if os.path.exists("hisat_mapping.finished"):
            pass
        else:
            self.hisat_mapping()
        if os.path.exists("sam2bam.finished"):
            pass
        else:
            self.sam_bam()
        if os.path.exists("transrate.finished"):
            pass
        else:
            self.run_transrate()
        self.set_output()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir = '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_denovo'
        data = {
            "id": "Transrate" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "denovo_rna_v2.transrate",
            "instant": True,
            "options": dict(
                assembly=test_dir + "/" + "Trinity.fasta",
                left=test_dir + "/" + "test_data1/reads.left.fq",
                right=test_dir + "/" + "test_data1/reads.right.fq",
                threads="8",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
