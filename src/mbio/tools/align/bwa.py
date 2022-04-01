# -*- coding: utf-8 -*-

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import glob
from biocluster.core.exceptions import OptionError


class BwaAgent(Agent):
    """
    bwa:比对工具
    version 1.0
    author: zhujuan
    last_modify: 2018.05
    """
    def __init__(self, parent):
        super(BwaAgent, self).__init__(parent)
        options = [
            {"name": "ref_database", "type": "string", "default": ""},  # 宿主参考序列库，eg:"plant，Brassica_napus"
            {"name": "ref_undefined", "type": "infile", "format": "sequence.fasta_dir"},
            # 未定义的宿主序列所在文件加，多个宿主cat到一个文件，并作为tool:align.bwa的输入文件
            {"name": "fq_type", "type": "string", "default": "PSE"},  # fq类型，PE、SE、PSE（即PE+SE，单端加双端）
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},
            # 输入质控后的fastq文件夹其中包含list文件
            {"name": "fastq_r", "type": "infile", "format": "sequence.fastq"},  # 右端序列文件
            {"name": "fastq_l", "type": "infile", "format": "sequence.fastq"},  # 左端序列文件
            {"name": "fastq_s", "type": "infile", "format": "sequence.fastq"},  # SE序列文件
            {"name": "head", "type": "string",
             "default": "'@RG\\tID:sample\\tLB:rna-seq\\tSM:sample\\tPL:ILLUMINA'"},  # 设置结果头文件
            {"name": "sam", "type": "outfile", "format": "align.bwa.sam_dir"},     # sam格式文件夹,内含对应list文件
            {"name": "method", "type": "string", "default": "align"},     # sam格式文件
            {"name": "result_path", "type": "string"},  # 指定bwa的生成文件存放路径，可不设置
            {"name": "db_path", "type": "string", "default": ""}, # 参考宿主库路径
        ]
        self.add_option(options)
        self._memory_increase_step = 30  # 每次重运行增加内存20G by gaohao 20190522

    def check_options(self):
        """
        检查参数
        """
        if self.option("ref_database") == "" and not self.option("ref_undefined").is_set:
            raise OptionError("请传入参考序列", code="31100701")
        if self.option("fastq_dir").is_set and not os.path.exists(self.option("fastq_dir").prop['path'] + "/list.txt"):
            raise OptionError("fastq序列文件夹需还有list文件", code="31100702")
        if self.option("method") == "align":
            if self.option('fq_type') not in ['PE', 'SE', 'PSE']:
                raise OptionError("请说明序列类型，PE or SE or 'PSE'?", code="31100703")
        if not self.option("fastq_dir").is_set and self.option('fq_type') in ["PE", "PSE"]:
            if not self.option("fastq_r").is_set:
                raise OptionError("请传入PE右端序列文件", code="31100704")
            if not self.option("fastq_l").is_set:
                raise OptionError("请传入PE左端序列文件", code="31100705")
        if not self.option("fastq_dir").is_set and self.option('fq_type') in ["SE", "PSE"]:
            if not self.option("fastq_s").is_set:
                raise OptionError("请传入SE序列文件", code="31100706")
        #if not self.option("fastq_dir").is_set:
        #    if self.option("result_path") == "":
        #        raise OptionError("请传入输出结果目录")
        return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 10
        # tmp_mem = 20 * (self._rerun_time + 1)  # 每次因拼接失败而重运行的内存增加20G by GHD @ 20180409
        # self._memory = '%sG' % tmp_mem
        self._memory = '20G'  # 改回 by guhaidong @ 20180427
        # self.logger.info('bwa use memory : ' + self._memory)


class BwaTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(BwaTool, self).__init__(config)
        self.bwa_path = "bioinfo/align/bwa-0.7.17/"
        self.ref_database = self.option("db_path") or self.config.SOFTWARE_DIR + "/database/Ens_Ref/"
        self.ref_fasta = ''
        if self.option("ref_database") != "":
            if os.path.exists(self.option("ref_database")):
                self.ref_fasta = self.option("ref_database")
            else:
                species_category = self.option("ref_database").replace(",", "/")
                self.ref_fasta = glob.glob(self.ref_database + species_category + "/*.sa")[0].rsplit(".sa")[0]
        if self.option("ref_undefined").is_set:
            ref_list = self.get_fastas()
            if len(ref_list) < 1:
                self.set_error("上传的宿主文件不存在！", code="31100701")
            elif len(ref_list) ==1:
                self.ref_fasta = ref_list[0]
            else:
                self.logger.info(ref_list)
                all_ref_undefined = os.path.join(self.option("ref_undefined").prop['path'],
                                                "ref_undefined/ref_undefined.fa")
                self.logger.info(all_ref_undefined)
                os.system('mkdir -p '+self.option("ref_undefined").prop['path'] + '/ref_undefined')
                os.system('cat ' + " ".join(ref_list) + ' >' + all_ref_undefined)
                self.ref_fasta = os.path.join(self.option("ref_undefined").prop['path'], "ref_undefined/ref_undefined.fa")
        if self.option("fastq_dir").is_set:
            self.samples = self.get_list()
            self.fq_dir = True
            self.fq_dir_path = self.option("fastq_dir").prop['path']
        else:
            self.fq_dir = False

    def get_fastas(self):
        all_fas = []
        for f in os.listdir(self.option("ref_undefined").prop["path"]):
            f_path = os.path.join(self.option("ref_undefined").prop["path"], f)
            f_split = f_path.split('.')
            if f_split[-1] in ['fa','fasta','fna']:
                all_fas.append(f_path)
            if f_split[-1] == 'gz' and f_split[-2] in ['fa','fasta','fna']:
                f_path2 = os.path.basename(f_path.rpartition('.gz')[0])
                os.system("rm {}*".format(f_path2))
                os.system("cp {} ./".format(f_path))
                cmd = 'gunzip ' + os.path.basename(f_path)
                if os.system(cmd) == 0:
                    all_fas.append(f_path2)
                else:
                    self.set_error("压缩的自定义宿主文件解压失败")
        return all_fas

    def bwa_index(self):
        cmd = "{}bwa index {}".format(self.bwa_path, self.ref_fasta)
        print cmd
        self.logger.info("开始构建参考序列索引")
        command = self.add_command("bwa_index", cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("成功构建参考序列索引！")
        elif command.return_code in [1, -7]:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180320
        else:
            self.logger.info("return code: %s" % command.return_code)
            self.set_error("构建索引出错", code="31100703")

    def bwa_aln(self, fastq, outfile):
        fq_name = fastq.split("/")[-1]
        cmd = "{}bwa aln -t 10 {} {} -f {}".format(self.bwa_path, self.ref_fasta, fastq, outfile)
        print(cmd)
        self.logger.info("开始运行{}_bwa_aln".format(fq_name.lower()))
        command = self.add_command("{}_bwa_aln".format(fq_name.lower()), cmd, ignore_error=True)
        command.run()
        return command

    def bwa_sampe(self, outfile, aln_l, aln_r, fastq_l, fastq_r):
        outfile_name = outfile.split("/")[-1]
        outfile = fastq_r.split("/")[-1].split(".")[0] + '.sam'
        if self.option("result_path") != "":
            outfile = os.path.join(self.option("result_path"), outfile)
        if self.option("head") == "":
            cmd = "{}bwa sampe -f {} {} {} {} {} {}".format(self.bwa_path, outfile,
                                                            self.ref_fasta, aln_l, aln_r, fastq_l, fastq_r)
        else:
            cmd = "{}bwa sampe -r {} -f {} {} {} {} {} {}".format(self.bwa_path, self.option("head"), outfile,
                                                                  self.ref_fasta, aln_l, aln_r, fastq_l, fastq_r)
        print(cmd)
        self.logger.info("开始运行{}_bwa_sampe".format(outfile_name.lower()))
        command = self.add_command("{}_bwa_sampe".format(outfile_name.lower()), cmd, ignore_error=True)
        command.run()
        if self.fq_dir is True:
            return command
        else:
            self.wait()
            if command.return_code == 0:
                self.logger.info("生成sam比对结果文件完成！")
            elif command.return_code == 1:
                self.logger.info("return code: %s" % command.return_code)
                self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180320
            else:
                self.logger.info("return code: %s" % command.return_code)
                self.set_error("生成sam比对结果文件出错", code="31100704")

    def bwa_samse(self, outfile, aln_s, fastq_s):
        outfile_name = outfile.split("/")[-1]
        outfile = fastq_s.split("/")[-1].split(".")[0] + '_s.sam'
        if self.option("result_path") != "":
            outfile = os.path.join(self.option("result_path"), outfile)
        if self.option("head") == "":
            cmd = "{}bwa samse -f {} {} {} {}".format(self.bwa_path, outfile, self.ref_fasta, aln_s, fastq_s)
        else:
            cmd = "{}bwa samse -r {} -f {} {} {} {}".format(self.bwa_path, self.option("head"),
                                                            outfile, self.ref_fasta, aln_s, fastq_s)
        print(cmd)
        self.logger.info("开始运行{}_bwa_sampe命令".format(outfile_name.lower()))
        command = self.add_command("{}_bwa_samse".format(outfile_name.lower()), cmd, ignore_error=True)
        command.run()
        if self.fq_dir is True:
            return command
        else:
            self.wait()
            if command.return_code == 0:
                self.logger.info("生成sam比对结果文件完成！")
            elif command.return_code == 1:
                self.logger.info("return code: %s" % command.return_code)
                self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180320
            else:
                self.logger.info("return code: %s" % command.return_code)
                self.set_error("生成sam比对结果文件出错", code="31100704")

    def multi_aln(self):
        samples = self.samples
        aln_commands = []
        for sample in samples:
            if self.option("fq_type") in ["PE", "PSE"]:
                aln_l_cmd = self.bwa_aln(os.path.join(self.fq_dir_path, samples[sample]["l"]),
                                         "{}_l.sai".format(sample))
                aln_r_cmd = self.bwa_aln(os.path.join(self.fq_dir_path, samples[sample]["r"]),
                                         "{}_r.sai".format(sample))
                aln_commands.append(aln_l_cmd)
                aln_commands.append(aln_r_cmd)
            if self.option("fq_type") in ["SE", "PSE"]:
                aln_s_cmd = self.bwa_aln(os.path.join(self.fq_dir_path, samples[sample]["s"]),
                                         "{}_s.sai".format(sample))
                aln_commands.append(aln_s_cmd)
                self.logger.info(aln_s_cmd)
        return aln_commands

    def multi_sam(self):
        samples = self.samples
        sam_commands = []
        sam_list_path = os.path.join(self.output_dir, "list.txt")
        with open(sam_list_path, "wb") as w:
            for sample in samples:
                if self.option("fq_type") in ["PE", "PSE"]:
                    fq_r = os.path.join(self.option("fastq_dir").prop['path'], samples[sample]["r"])
                    fq_l = os.path.join(self.option("fastq_dir").prop['path'], samples[sample]["l"])
                    sam_pe_cmd = self.bwa_sampe(os.path.join(self.output_dir, "{}.sam".format(sample)),
                                                "{}_l.sai".format(sample), "{}_r.sai".format(sample), fq_l, fq_r)
                    sam_commands.append(sam_pe_cmd)
                    w.write(sample+".sam\t"+sample+"\tpe\n")
                if self.option("fq_type") in ["SE", "PSE"]:
                    fq_s = os.path.join(self.option("fastq_dir").prop['path'], samples[sample]["s"])
                    sam_se_cmd = self.bwa_samse(os.path.join(self.output_dir, "{}_s.sam".format(sample)),
                                                "{}_s.sai".format(sample), fq_s)
                    sam_commands.append(sam_se_cmd)
                    w.write(sample+"_s.sam\t"+sample+"\tse\n")
            return sam_commands

    def get_list(self):
        list_path = self.option("fastq_dir").prop['path'] + "/list.txt"
        self.logger.info(list_path)
        list_path = os.path.join(self.option("fastq_dir").prop['path'], "list.txt")
        if os.path.exists(list_path):
            self.logger.info(list_path)
        sample = {}
        with open(list_path, "rb") as l:
            for line in l:
                line = line.strip().split()
                if len(line) == 3:
                    if line[1] not in sample:
                        sample[line[1]] = {line[2]: line[0]}
                    else:
                        sample[line[1]][line[2]] = line[0]
                if len(line) == 2:
                    if line[1] not in sample:
                        sample[line[1]] = {"s": line[0]}
                    else:
                        sample[line[1]]["s"] = line[0]
        return sample

    def set_ouput(self):
        self.logger.info('开始设置输出结果文件')
        try:
            self.option('sam', self.output_dir)
            self.logger.info("设置输出结果文件正常")
        except Exception as e:
            self.set_error("设置输出结果文件异常——%s", variables=(e), code="31100705")

    def bwa_mem(self, read, read2=None):
        fq_name = read.split('/')[-1].split('.')[0]
        outfile = fq_name
        self.logger.info('### {} ###'.format(self.option("result_path")))
        if self.option("result_path"):
            self.logger.info('### {} ###'.format(self.option("result_path")))
            outfile = os.path.join(self.option("result_path"), outfile)
        cmd = os.path.join(self.config.SOFTWARE_DIR, self.bwa_path) + 'bwa mem -t 10 '
        if self.option('head'):
            cmd += '-R {} '.format(self.option('head'))
        cmd += '{} {}'.format(self.ref_fasta, read)
        if read2:
            cmd += ' {} > {}.sam'.format(read2, outfile)
        else:
            cmd += ' > {}_s.sam'.format(outfile)
        command = self.add_command("{}_bwa_mem".format(fq_name.lower()), cmd, shell=True, ignore_error=True).run()
        return command

    def multi_mem(self):
        samples = self.samples
        mem_commands = []
        for sample in self.samples:
            if self.option("fq_type") in ["PE", "PSE"]:
                pe_mem = self.bwa_mem(self.option("fastq_l").prop['path'], self.option("fastq_r").prop['path'])
                mem_commands.append(pe_mem)
            if self.option("fq_type") in ['SE', 'PSE']:
                s_mem = self.bwa_mem(self.option("fastq_s").prop['path'])
                mem_commands.append(s_mem)
        return mem_commands

    def run(self):
        super(BwaTool, self).run()
        if self.option('method') == 'index':
            self.bwa_index()
        else:
            if os.path.exists(self.ref_fasta + '.sa'):
                pass
            else:
                self.bwa_index()
            if self.option('fastq_dir').is_set:
                mem_commands = self.multi_mem()
                self.logger.info(mem_commands)
                self.wait()
                for mem_cmd in mem_commands:
                    if mem_cmd.return_code == 0:
                        self.logger.info('运行{}完成'.format(mem_cmd.name))
                    elif mem_cmd.return_code in [137, 135, 1]:
                        self.logger.info("return code: %s" % mem_cmd.return_code)
                        self.add_state('memory_limit', 'memory is low!')
                    else:
                        self.logger.info("return code: %s" % mem_cmd.return_code)
                        self.set_error("运行%s运行出错!", variables=(mem_cmd.name), code="31100706")
                        return False
            else:
                if self.option("fq_type") in ["PE", "PSE"]:
                    pe_mem = self.bwa_mem(self.option("fastq_l").prop['path'], self.option("fastq_r").prop['path'])
                    self.wait(pe_mem)
                    if pe_mem.return_code == 0:
                        self.logger.info('pe_mem 比对完成')
                    elif pe_mem.return_code in [137, 135, 1]:
                        self.logger.info("return code: %s" % pe_mem.return_code)
                        self.add_state('memory_limit', 'memory is low!')
                    else:
                        self.set_error('比对出错')
                if self.option("fq_type") in ['SE', 'PSE']:
                    s_mem = self.bwa_mem(self.option("fastq_s").prop['path'])
                    self.wait(s_mem)
                    if s_mem.return_code == 0:
                        self.logger.info('s_mem 比对完成')
                    elif s_mem.return_code in [137, 135, 1]:
                        self.logger.info("return code: %s" % s_mem.return_code)
                        self.add_state('memory_limit', 'memory is low!')
                    else:
                        self.set_error('比对出错')
        self.set_ouput()
        self.end()

    '''
    def run(self):
        """
        运行
        """
        super(BwaTool, self).run()
        if self.option("method") == "index":
            self.bwa_index()
        else:
            if os.path.exists(self.ref_fasta + ".sa"):
                pass
            else:
                self.bwa_index()
            if self.option("fastq_dir").is_set:
                aln_commands = self.multi_aln()
                self.logger.info(aln_commands)
                self.wait()
                for aln_cmd in aln_commands:
                    if aln_cmd.return_code == 0:
                        self.logger.info("运行{}完成".format(aln_cmd.name))
                    elif aln_cmd.return_code == 1:
                        self.logger.info("return code: %s" % aln_cmd.return_code)
                        self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180320
                    else:
                        self.logger.info("return code: %s" % aln_cmd.return_code)
                        self.set_error("运行%s运行出错!", variables=(aln_cmd.name), code="31100706")
                        return False
                sam_commands = self.multi_sam()
                self.logger.info(sam_commands)
                self.wait()
                for sam_cmd in sam_commands:
                    if sam_cmd.return_code == 0:
                        self.logger.info("运行{}完成".format(sam_cmd.name))
                    elif sam_cmd.return_code == 1:
                        self.logger.info("return code: %s" % sam_cmd.return_code)
                        self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180320
                    else:
                        self.logger.info("return code: %s" % sam_cmd.return_code)
                        self.set_error("运行%s运行出错!", variables=(sam_cmd.name), code="31100707")
                        return False
            else:
                if self.option("fq_type") in ["PE", "PSE"]:
                    aln_l = self.bwa_aln(self.option("fastq_l").prop['path'], "aln_l.sai")
                    aln_r = self.bwa_aln(self.option("fastq_r").prop['path'], "aln_r.sai")
                    self.wait(aln_l, aln_r)
                    if aln_l.return_code == 0:
                        self.logger.info("左端比对完成！")
                    elif aln_l.return_code in [1, -9]:  # bwa会因内存不足导致比对出错 modified by GHD @ 20180515 @ 20190128
                        self.logger.info("return code: %s" % aln_l.return_code)
                        self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180320
                    else:
                        self.logger.info("return code: %s" % aln_l.return_code)
                        self.set_error("左端比对出错", code="31100708")
                    if aln_r.return_code == 0:
                        self.logger.info("右端比对完成！")
                    elif aln_r.return_code in [1, -9]:  # bwa会因内存不足导致比对出错 modified by GHD @ 20180515 @ 20190128
                        self.logger.info("return code: %s" % aln_r.return_code)
                        self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180320
                    else:
                        self.logger.info("return code: %s" % aln_r.return_code)
                        self.set_error("右端比对出错", code="31100709")
                    self.bwa_sampe("pe.sam", "aln_l.sai", "aln_r.sai", self.option("fastq_l").prop['path'],
                                   self.option("fastq_r").prop['path'])
                elif self.option("fq_type") in ["SE", "PSE"]:
                    aln_s = self.bwa_aln(self.option("fastq_s").prop['path'], "aln_s.sai")
                    self.wait(aln_s)
                    if aln_s.return_code == 0:
                        self.logger.info("比对完成！")
                    elif aln_s.return_code in [1, -9]:  # modified by ghd @ 20190128
                        self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20190128
                    else:
                        self.set_error("比对出错", code="31100710")
                    self.bwa_samse("se.sam", "aln_s.sai", self.option("fastq_s").prop['path'])
        self.set_ouput()
        self.end()
        '''
