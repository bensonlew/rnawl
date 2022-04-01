# -*- coding: utf-8 -*-
# __author__ = 'guoquan'  
from .admin import Admin
import sys
import importlib
import os
import re
reload(sys)
sys.setdefaultencoding('utf8')


class PymoduleAction(Admin):

    def __init__(self):
        super(PymoduleAction, self).__init__()
        self.nav_index = 0
        self.title = "Python代码显示"
        self.check_required(["name", "path"], relation="or")
        self.file_path = None
        self._module_path = None

    def GET(self):
        if self.params.name and self.params.name != "":
            try:
                importlib.import_module(self.params.name)
                if self.params.name in sys.modules.keys():
                    path = sys.modules[self.params.name].__file__
                    if re.search(r"pyc$", path):
                        self.file_path = re.sub(r"pyc$", "py", path)
                    else:
                        self.file_path = path
                self._module_path = self.params.name
            except:
                return self.failed("模块无法加载或者不存在！")
        else:
            self.file_path = os.path.join(os.environ['HOME'], self.params.path)

        if not os.path.exists(self.file_path):
            return self.failed("模块无法加载！")
        return self.render.admin.pymodule(self)

    @property
    def code(self):
        data = ""
        with open(self.file_path, "r") as f:
            line = f.readline()
            while line:
                m = re.match(r"^\s*from\s+(\S+)(.*)$", line)
                if m:
                    data += "from <a target=\"_blank\" href=\"%s\">%s</a> %s\n" % (self.get_module_path(m.group(1)),
                                                                         m.group(1), m.group(2))
                else:
                    mm = re.match(r"^\s*import\s+(\S+)\s*$", line)
                    if mm:
                        data += "import <a target=\"_blank\" href=\"%s\">%s</a>\n" % (
                            self.get_module_path(mm.group(1)), mm.group(1))
                    else:
                        data += line
                line = f.readline()
        return data

    def get_module_path(self, module):

        m1 = re.match(r"^(\.+)(.+)$", module)
        if m1:
            dots = m1.group(1)
            rel_path = m1.group(2)
            if self._module_path:
                path_list = self._module_path.split(".")
                if len(path_list) < len(dots):
                    return "pymodule?name=%s" % module
                else:
                    for i in range(len(dots)):
                        path_list.pop(-1)
                    path_list.extend(rel_path.split("."))
                return "pymodule?name=%s" % ".".join(path_list)
            else:
                path_list = self.file_path.split("/")
                if len(path_list) < len(dots):
                    return "pymodule?name=%s" % module
                else:
                    for i in range(len(dots)):
                        path_list.pop(-1)
                    path_list.extend(rel_path.split("."))
                return "pymodule?path=%s.py" % "/".join(path_list)
        else:
            return "pymodule?name=%s" % module

    @property
    def line(self):
        if self.params.line:
            try:
                line = int(self.params.line)
            except:
                line = 0
            return line
        else:
            return 0
