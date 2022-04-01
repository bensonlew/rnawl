# -*- coding: utf-8 -*-
# __author__ = 'guoquan'  
import web
import urlparse
import urllib
import os


class Paginatiton(object):

    def __init__(self, model):
        self._model = model
        self._url = urlparse.urlparse(web.ctx.fullpath)
        self._current_page = self._model.current_page

    @property
    def current_page(self):
        return self._current_page

    def url(self, page):
        new_query = []
        has_query = False
        if self._url.query:
            query = urlparse.parse_qsl(self._url.query)
            for q in query:
                if q[0] == "page":
                    has_query = True
                    new_query.append(("page", page))
                else:
                    new_query.append(q)
        if not has_query:
            new_query.append(("page", page))
        query_str = urllib.urlencode(new_query)
        new_parse_result = urlparse.ParseResult(scheme=self._url.scheme, netloc=self._url.netloc, path=self._url.path,
                                                params=self._url.params, query=query_str, fragment=self._url.fragment)
        return urlparse.urlunparse(new_parse_result)

    def get_pager(self):
        path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../views/'))
        render = web.template.render(path)
        return render.admin.pagination(self)

    @property
    def page_list(self):
        if self._model.total_page < 11:
            return [i for i in range(1, self._model.total_page + 1)]
        elif self.current_page <= self._model.total_page - 5:
            if self.current_page > 5:
                return [i for i in range(self.current_page - 4, self.current_page + 6)]
            else:
                return [i for i in range(1, 11)]
        else:
            return [i for i in range(self._model.total_page - 9, self._model.total_page + 1)]

    @property
    def last(self):
        return self._model.total_page

    @property
    def previous(self):
        if self._current_page > 1:
            return self._current_page - 1
        else:
            return 1

    @property
    def next_page(self):
        if self._current_page < self.last:
            return self._current_page + 1
        else:
            return self.last
