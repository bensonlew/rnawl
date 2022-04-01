# -*- coding: utf-8 -*-
# __author__ = 'guoquan'  
from .admin import Admin


class ApiAction(Admin):
    def __init__(self):
        super(ApiAction, self).__init__()
        self.nav_index = 9
        self.title = "API日志"
        self.list_data = []
        self._table = "apilog"

    def GET(self):
        self.list_model.pagesize = 50
        where_str = None

        if self.params.data_type:
            if self.params.table_search:
                if self.params.data_type == 'task_id':
                    where_str = "task_id='%s'" % self.params.table_search
                if self.params.data_type == 'api':
                    where_str = "api='%s'" % self.params.table_search
            elif self.params.time_from or self.params.time_to:
                key = "addtime"
                if self.params.time_from:
                    where_str = "%s>='%s' " % (key, self.params.time_from)
                if self.params.time_to:
                    if where_str:
                        where_str += "and %s<='%s' " % (key, self.params.time_to)
                    else:
                        where_str = "%s<='%s' " % (key, self.params.time_to)
        self.list_model.set_conditions(where=where_str)
        if self.list_model.total > 0:
            for d in self.list_model.get_data(self.current_page):
                self.list_data.append(d)
        return self.render.admin.api(self)

