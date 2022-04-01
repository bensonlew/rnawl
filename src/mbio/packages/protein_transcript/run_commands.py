# -*- coding: utf-8 -*-
# author fengyitong 2019-01

import os
import logging
import random
import json
import time


class RunCommands(object):
    # def __init__(self, cmd = None, name = None, run_wd = os.getcwd(), qsub = True, memory = 1, node = 1):
    def __init__(self, params = None):
        LOG_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
        self.log = logging
        if not 'name' in params or not params['name']:
            self.name = 'command' + str(random.randint(0,10000))
        else:
            self.name = params['name']
        self.cmd = None
        if 'cmd' in params:
            self.cmd = params['cmd']
        self.qsub = 1
        if 'qsub' in params:
            self.qsub = params['qsub']
        self.memory = 1
        if 'memory' in params:
            self.memory = params['memory'].strip('G')
        self.node = 1
        if 'node' in params:
            self.node = params['node']
        run_wd = os.getcwd()
        if 'run_wd' in params:
            run_wd = params['run_wd']
#----judge to qsub or slurm----
        if u'/centos7users/' in run_wd and not 'qsub' in params:
            self.qsub = 2
        self.group = 'protein'
        try:
            group = os.popen('hostname').read().split('.')[0]
        except:
            pass
        if 'group' in params:
            self.group = params['group']
        self.work_dir = os.path.join(run_wd, self.name)
        if not os.path.exists(self.work_dir):
            os.mkdir(self.work_dir)
        self.end = 0
        self.level = 0
        # self.log.basicConfig(filename=os.path.join(run_wd, self.name + '.log'), level=logging.INFO, format=LOG_FORMAT)
        self.logfile = os.path.join(run_wd, 'all_run.log')
        self.log.basicConfig(filename=self.logfile, level=logging.DEBUG, format=LOG_FORMAT)
        self.running = os.path.join(self.work_dir, '.' + self.name + '.running')
        self.judge = os.path.join(self.work_dir, '.' + self.name + '.finished')

    def set_params(self, params):
        if not type(params) == dict or not 'cmd' in params:
            self.end(normal=1,out='your command is not standardily %s'%json.dumps(params))
        if 'cmd' in params:
            self.cmd = params['cmd']
        if 'qsub' in params:
            self.qsub = params['qsub']
        if 'memory' in params:
            self.memory = params['memory']
        if 'node' in params:
            self.node = params['node']
        if 'group' in params:
            self.group = params['group']

    def run(self):
        if os.path.exists(self.running):
            return
        if os.path.exists(self.judge):
            return
        self.log.info("The command %s start running, the command info is %s"%(self.name, self.cmd))
        # os.system('touch ' + self.running)
        error_path = os.path.join(self.work_dir, self.name + '.e')
        out_path = os.path.join(self.work_dir, self.name + '.o')
        if self.qsub == 1:
            os.system('touch ' + self.running)
            cmd_str = '#PBS -l nodes=1:ppn=%s\n#PBS -l mem=%s\n#PBS -q %s\n'%(str(self.node), str(self.memory)+'G', self.group)
            cmd_str += '#PBS -o %s\n#PBS -e %s\n'%(out_path, error_path)
            cmd_str += 'cd ' + self.work_dir + '\n'
            cmd_str += self.cmd + '\n'
            cmd_str += '/bin/rm ' + self.running + '\n'
            cmd_str += 'touch ' + self.judge + '\n'
            bash_ = os.path.join(self.work_dir, 'qsub_' + self.name + '.bash')
            with open(bash_, 'w') as fw:
                fw.write(cmd_str)
            qsub_id = os.popen('qsub ' + bash_).read().split('.')[0]
            # error_path = "unknown"
            # out_path = "unknown"
            # for line in os.popen('qstat -f ' + qsub_id).readlines():
            #     if u'Error_Path = ' in line:
            #         error_path = line.strip().split(':')[-1]
            #     if u'Output_Path = ' in line:
            #         out_path = line.strip().split(':')[-1]
            self.log.debug("The command %s has been qsub, the error_path is %s, the out_path is %s" % (
                self.name, error_path, out_path))
        elif self.qsub == 2:
            os.system('touch ' + self.running)
            node = os.popen('hostname').read().split('.')[0]
            cmd_str = '#!/bin/bash\n#SBATCH -D %s\n#SBATCH -N 1\n#SBATCH --ntasks-per-node %s\n#SBATCH -J %s\n#SBATCH -p %s\n#SBATCH --mem=%sG\n#SBATCH -o %s\n#SBATCH -e %s\n' % (str(self.work_dir), str(self.node), self.name, self.group, self.memory, out_path, error_path)
            cmd_str += 'cd ' + self.work_dir + '\n'
            cmd_str += self.cmd + '\n'
            cmd_str += '/bin/rm ' + self.running + '\n'
            cmd_str += 'touch ' + self.judge + '\n'
            bash_ = os.path.join(self.work_dir, 'sbatch_' + self.name + '.sbatch')
            with open(bash_, 'w') as fw:
                fw.write(cmd_str)
            if node != self.group:
                tmp = os.getcwd()
                c_tmp = 'ssh -Y ' + node
                c_tmp_ = 'cd ' + tmp
                os.system(c_tmp)
                os.system(c_tmp_)
            os.system('sbatch ' + bash_)
            self.log.debug("The command %s has been qsub, the error_path is %s, the out_path is %s" % (
                self.name, error_path, out_path))
        elif self.qsub == 0:
            os.system('touch ' + self.running)
            cmd_str = '#!/bin/bash\n'
            cmd_str += 'cd ' + self.work_dir + '\n'
            cmd_str += self.cmd + '\n'
            cmd_str += '/bin/rm ' + self.running + '\n'
            cmd_str += 'touch ' + self.judge + '\n'
            bash_ = os.path.join(self.work_dir, 'nohup_' + self.name + '.sh')
            with open(bash_, 'w') as fw:
                fw.write(cmd_str)
            os.system('nohup bash ' + bash_ + ' > %s 2>&1 &' % out_path)
            self.log.debug("The command %s has been qsub, the error_path is %s, the out_path is %s" % (
                self.name, error_path, out_path))
        else:
            # 在centos7上无法slurm的对象再调用slurm，只能本地跑，本地跑的话，如果两个函数一起切换目录就会出问题，所以给他睡一会
            time.sleep(random.random())
            self.log.debug("The command %s will been runned on local site\n first, we change dir to %s" % (
            self.name, self.work_dir))
            tmp = os.getcwd()
            os.chdir(self.work_dir)
            time.sleep(2)
            cmds = self.cmd.split('\n')
            while '' in cmds:
                cmds.remove('')
            for cmd in cmds:
                with open(self.logfile, 'r') as lr:
                    if '_____'+cmd+'______' in lr.read():
                        print(cmd+' had runned')
                        continue
                cmd_info = os.popen(cmd).read()
                self.log.debug('The command %s now is running _____%s______\n the output in %s' %(self.name, cmd, cmd_info))
            else:
                # os.system('/bin/rm ' + self.running)
                os.system('touch ' + self.judge)
                os.chdir(tmp)
                self.log.debug('The command %s finished'%self.name)

    def end(self, normal=0, out=None):
        if not normal:
            print('your program finished normally')
        else:
            print('your program finished unexceptally')
        if out:
            print(out)
        os._exit(normal)