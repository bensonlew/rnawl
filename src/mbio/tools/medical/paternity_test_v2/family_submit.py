# !/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python
# -*- coding: utf-8 -*-
# __auther__ = 'hongdong'
from biocluster.core.exceptions import OptionError
from socket import error as SocketError
from biocluster.agent import Agent
from biocluster.tool import Tool
import httplib
import urllib2
import urllib
import hashlib
import random
import time
import json
import copy


class FamilySubmitAgent(Agent):
    """
    该tool用于将组建好的不同家系按照接口形式投递出去运行，输入是一个家系二维列表[[WQ12345678-F,WQ12345678-M,WQ12345678-S]]
    version v1.0
    author: hongdongxuan
    last_modify: 20171117
    """
    def __init__(self, parent):
        super(FamilySubmitAgent, self).__init__(parent)
        options = [
            {"name": "family_list", "type": "string"},  # 输入的是家系列表
            {"name": "datasplit_id", "type": "string"},  # 拆分批次id，一个拆分表就是一个批次
            {"name": "types", "type": "string", "default": 'sanger'},  # 测试环境类型，tsg，tsanger，sanger
            {"name": "member_id", "type": "string"},
            {"name": "fastq_path", "type": "string"},
            {"name": "err_min_num", "type": "int", "default": 9}
        ]
        self.add_option(options)
        self.step.add_steps("family_submit")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)
        self._cpu = None
        self._memory = None

    def stepstart(self):
        self.step.family_submit.start()
        self.step.update()

    def stepfinish(self):
        self.step.family_submit.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("family_list"):
            raise OptionError("必须输入家系列表")
        if not self.option("fastq_path"):
            raise OptionError("必须输入fastq文件所在路径")
        if not self.option("member_id"):
            raise OptionError("必须输入会员id")
        if not self.option("datasplit_id"):
            raise OptionError('必须输入拆分的批次id')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(FamilySubmitAgent, self).end()


class FamilySubmitTool(Tool):
    """
    运行获取客户信息，然后导入到mongo表中
    """
    def __init__(self, config):
        super(FamilySubmitTool, self).__init__(config)
        self._version = '1.0.1'
        self.web_data = None
        self.url = None

    def webapi_submit_run(self):
        """
        循环投递出不同的家系分析的接口，逻辑待写
        """
        self.web_data = json.loads(self.option("family_list"))
        for family in self.web_data:
            # if not self.api.api("medical.paternity_test_v2").tab_exist(family[0]):
            #     raise Exception("父本{}tab文件不存在！".format(family[0]))
            # else:
            #     self.logger.info("父本{}tab文件存在！".format(family[0]))
            # if not self.api.api("medical.paternity_test_v2").tab_exist(family[1]):
            #     raise Exception("母本{}tab文件不存在！".format(family[1]))
            # else:
            #     self.logger.info("母本{}tab文件存在！".format(family[1]))
            params = {"dad_id": family[0], "mom_id": family[1], "son_id": family[2]}
            Submit(self, params, self.option("types")).post_to_webapi()
            time.sleep(10)

    def run(self):
        super(FamilySubmitTool, self).run()
        self.webapi_submit_run()
        self.end()


class Submit(object):
    """设定投递对象"""
    def __init__(self, bind_object, params, types):
        if types == 'tsg':
            self._url = "http://10.101.203.193:9090"  # 这里用内网网段，这样投递到集群中才能获取到mysql，192的那个不行
            self._client = "client03"
            self.mysql_client_key = 'hM4uZcGs9d'
        elif types == 'tsanger':
            self._url = "http://bcl.tsanger.com"
            self._client = "client03"
            self.mysql_client_key = "hM4uZcGs9d"
        elif types == 'sanger':
            self._url = "http://bcl.i-sanger.com"
            self._client = "client01"
            self.mysql_client_key = "1ZYw71APsQ"
        else:
            raise Exception("测试机器的类型{}不正确！".format(types))
        self._api = "/s/pt/family_analysis"
        self._params = params
        self.bind_object = bind_object
        self.json_params = None

    def signature(self):
        timestamp = str(int(time.time()))
        nonce = str(random.randint(1000, 10000))
        web_key = self.mysql_client_key
        sha1 = hashlib.sha1()
        key_list = [web_key, nonce, timestamp]
        key_list.sort()
        map(sha1.update, key_list)
        hashkey = sha1.hexdigest()
        signature = {
            "client": self._client,
            "nonce": nonce,
            "timestamp": timestamp,
            "signature": hashkey
        }
        return urllib.urlencode(signature)

    def url_params_format(self):
        temp_params = copy.deepcopy(self._params)
        temp_params['batch_id'] = str(self.bind_object.option('datasplit_id'))
        temp_params['fastq_path'] = str(self.bind_object.option("fastq_path"))
        temp_params['err_min_num'] = str(self.bind_object.option("err_min_num"))
        temp_params['member_id'] = str(self.bind_object.option("member_id"))
        return temp_params

    def post(self):
        """用于投递接口"""
        httpHandler = urllib2.HTTPHandler(debuglevel=1)
        httpsHandler = urllib2.HTTPSHandler(debuglevel=1)
        opener = urllib2.build_opener(httpHandler, httpsHandler)
        urllib2.install_opener(opener)
        params = self.url_params_format()
        dataa = urllib.urlencode(params)
        url = "%s?%s" % (self._url + self._api, self.signature())
        request = urllib2.Request(url, dataa)
        for i in range(3):
            try:
                response = urllib2.urlopen(request)
            except (urllib2.HTTPError, urllib2.URLError, httplib.HTTPException, SocketError) as e:
                self.bind_object.logger.warn("接口投递失败, 第{}次(最多三次尝试), 错误:{}".format(i + 1, e))
            else:
                the_page = response.read()
                return json.loads(the_page)
        return {"info": "拆分接口投递失败！", "success": False}

    def params_pack(self, dict_params):
        """参数打包"""
        return json.dumps(dict_params, sort_keys=True, separators=(',', ':'))

    def check_end(self):
        """检查家系分析投递出去运行的状态，暂时用不到"""
        pass

    def check_params(self):
        """检查参数一致，检车然后确定是否要投递接口，参数一致就不投递了！,
        返回值为true说明找不到，会进行投递分析，否则就pass"""
        self.json_params = self.params_pack(self._params)
        return self.bind_object.api.api("medical.paternity_test_v2").find_params_result(self.json_params)

    def post_to_webapi(self):
        """投递接口"""
        if self.check_params():
            result = self.post()
            if result['success']:
                self.bind_object.logger.info("params:{}投递成功！{}".format(self._params, result['content']))
            else:
                self.bind_object.logger.info("params:{}投递失败！{}".format(self._params, result['info']))
        else:
            self.bind_object.logger.info("参数已经存在不进行投递接口分析！")
            # self.bind_object.api.api("medical.paternity_test_v2").update_process(self.bind_object.option("datasplit_id"))
