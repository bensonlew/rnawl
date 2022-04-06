# !/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/python
# -*- coding: utf-8 -*-
# __auther__ = 'yuguo'


from __future__ import print_function
from biocluster.core.exceptions import OptionError
from socket import error as SocketError
from biocluster.agent import Agent
from biocluster.tool import Tool
import httplib
import gevent
import urllib2
import urllib
import hashlib
import random
import time
import json


class SubmitWorkAgent(Agent):
    """
    该tool用于管理医学样本报告分析投递顺序。
    注： failed_to_stop参数，针对流程不管投递失败还是计算失败流程停止， 如果为true的时候，不管投递的对象的运行状态，
    流程都不会停止，默认出错是不会停止的
    """
    def __init__(self, parent):
        super(SubmitWorkAgent, self).__init__(parent)
        options = [
            {"name": "api_name", "type": "string"},
            {"name": "params", "type": "string"},
            {"name": "failed_to_stop", "type": "string", "default": "false"},   # 判断该tool运行失败是否终止整个流程
            {"name": "is_keep_process", "type": "string", "default": "true"}  # 用于处理是不是要实时监控tool的运行状态
            # ，还是直接投递成功就ok了
        ]
        self.add_option(options)
        self.queue = "gada"

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("api_name"):
            raise OptionError("缺少必要参数api_name：接口名称")
        if not self.option("params"):
            raise OptionError("缺少必要参数params：接口参数")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = '1G'

    def end(self):
        super(SubmitWorkAgent, self).end()


class SubmitWorkTool(Tool):
    """
    """
    def __init__(self, config):
        super(SubmitWorkTool, self).__init__(config)
        self.workinfo = ""

    def submit_task(self):
        self.logger.info("开始提交，api_name:{},params:{}".format(self.option('api_name'), self.option('params')))
        task = Submit(self, self.option('api_name'), json.loads(self.option('params')), 'tsanger')
        result = task.post_to_webapi()
        if self.option("is_keep_process") == "true":
            if result.split(":")[0] == 'Failed':
                self.workinfo = result.split(":")[0]
                if self.option("failed_to_stop") == 'true':
                    self.set_error('Tool:submit_work运行出错，{}'.format(result))
                    raise Exception("Tool:submit_work运行出错，{}".format(result))
                else:
                    self.logger.info('Tool:submit_work运行出错：{}'.format(result))
            else:
                self.workinfo = result.split("::")[1]
                if self.workinfo.split(':')[0] == 'failed':
                    if self.option("failed_to_stop") == 'true':
                        self.set_error('Tool:submit_work正常结束，投递的workflow运行失败：{}'.format(self.workinfo))
                        raise Exception('Tool:submit_work正常结束，投递的workflow运行失败：{}'.format(self.workinfo))
                    else:
                        self.logger.info('Tool:submit_work正常结束，投递的workflow运行失败：{}'.format(self.workinfo))
                else:
                    self.logger.info('Tool:submit_work正常结束，投递的workflow运行成功：{}'.format(self.workinfo))
        else:
            self.logger.info(result)

    def run(self):
        super(SubmitWorkTool, self).run()
        self.submit_task()
        # self.remote.sync(self.workinfo)
        self.end()


class Submit(object):
    def __init__(self, bind_object, api_name, params, types):
        """
        __auther__ = 'hongdong
        用于投递接口  modified at 20171130,
        :param bind_object:
        :param api_name: 接口投递的url，现在是自动加载的
        :param params: 接口需要的参数
        :param types:tsg, tsanger, sanger,后面待改进
        """
        if types == 'tsg':
            self._url = "http://bcl.tsg.com"  # 这里用内网网段，这样投递到集群中才能获取到mysql，192的那个不行
            # self._url = "http://192.168.12.102:9090"  # 这里用内网网段，这样投递到集群中才能获取到mysql，192的那个不行
            self._client = "client03"
            self.mysql_client_key = 'hM4uZcGs9d'
        elif types == 'tsanger':
            self._url = "http://bcl.tsanger.com"
            self._client = "client03"
            self.mysql_client_key = "hM4uZcGs9d"
        elif types == 'sanger':
            self._url = "http://bcl.sanger.com"
            self._client = "client01"
            self.mysql_client_key = "1ZYw71APsQ"
        else:
            raise Exception("测试机器的类型{}不正确！".format(types))
        self._api = api_name
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

    def post(self):
        """
        用于投递接口
        :return:
        """
        httpHandler = urllib2.HTTPHandler(debuglevel=1)
        httpsHandler = urllib2.HTTPSHandler(debuglevel=1)
        opener = urllib2.build_opener(httpHandler, httpsHandler)
        urllib2.install_opener(opener)
        # params = self.set_params()
        dataa = urllib.urlencode(self._params)
        url = "%s/%s?%s" % (self._url, self._api, self.signature())
        request = urllib2.Request(url, dataa)
        for i in range(3):
            try:
                response = urllib2.urlopen(request)
            except (urllib2.HTTPError, urllib2.URLError, httplib.HTTPException, SocketError) as e:
                gevent.sleep(60)
                self.bind_object.logger.warn("接口投递失败, 第{}次(最多三次尝试), 错误:{}".format(i + 1, e))
            else:
                the_page = response.read()
                return json.loads(the_page)
        return {"info": "拆分接口投递失败！", "success": False}

    def check_status(self, id_, collection):
        """
        用于检查该分析是否完成, 分析状态分为start，end，failed
        :param id_:
        :param collection:
        :return:
        """
        stat = ""
        while True:
            status = self.bind_object.api.api("medical.paternity_test_v3.paternity_test_v3").check_end(id_, collection)
            if str(status) in ["end", "failed"]:
                # self.bind_object.logger.info("接口对应workflow运行成功！")
                stat = str(status)
                break
            else:
                gevent.sleep(60)
        return stat

    def post_to_webapi(self):
        """
        投递接口,接口投递状态分为投递成功与投递失败
        :return:
        """
        self.bind_object.logger.info("开始post_to_webapi")
        result = self.post()
        stat = True
        if result['success']:
            self.bind_object.logger.info("投递接口成功！params:{},{},wo"
                                         "rkflow_id({})".format(self._params, result['content'], result['workflow_id']))
            if self.bind_object.option("is_keep_process") == "true":
                stat = self.check_status(result['content']['ids']['id'], result['content']['ids']['name'])
            return 'Success:投递接口成功::{}:workflow_id({})'.format(stat, result['workflow_id'])
        else:
            self.bind_object.logger.info("投递接口失败！params:{},{}".format(self._params, result['info']))
            return 'Failed:投递接口失败:{}'.format(result['info'])
