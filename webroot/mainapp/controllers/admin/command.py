# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
from .admin import Admin
import re
from biocluster.core.function import friendly_size


class CommandAction(Admin):
    def __init__(self):
        super(CommandAction, self).__init__()
        self.nav_index = 8
        self.title = "Command命令"
        self.list_data = []

    def GET(self):
        self.list_model.pagesize = 50
        if self.list_model.total > 0:
            for d in self.list_model.get_data(self.current_page):
                setattr(d, "cpu_used", d.max_cpu_use / 100)
                setattr(d, "mem_used", friendly_size(d.max_vms))
                setattr(d, "avg_cpu_used", d.average_cpu_use / 100)
                setattr(d, "avg_mem_used", friendly_size(d.average_vms))
                self.list_data.append(d)
        return self.render.admin.command(self)

    def memory(self, value):
        """
        返回Job申请的内存大小，单位为字节数
        :return:
        """
        if value == "":
            return 0
        pattern = re.compile(r'([\d\.]+)([gmk])', re.I)
        match = re.match(pattern, value)
        if match:
            unit = match.group(2)
            if unit.upper() == "G":
                return float(match.group(1)) * 1024 * 1024 * 1024
            elif unit.upper() == "M":
                return float(match.group(1)) * 1024 * 1024
            elif unit.upper() == "K":
                return float(match.group(1)) * 1024
        else:
            pattern = re.compile(r'([\d\.]+)', re.I)
            match = re.match(pattern, value)
            if match:
                return float(match.group(1))
        return 0
