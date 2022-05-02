# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
from .ssh import SSH
import random


class SSH1(SSH):
    """
    扩展通过配置文件中指定的IP投递任务
    """
    def __init__(self, agent):
        super(SSH1, self).__init__(agent)
        self.config = agent.config
        self.server_ip = self.get_remote_ip()

    def get_remote_ip(self):
        """
        获取远程主机IP

        :return:
        """
        mode = self.config.SSH1_MODE
        ip_list = self.config.SSH1_IP_LIST
        if mode == "random":
            return random.choice(ip_list)

