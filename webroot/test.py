# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

import web
import mainapp.core.auto_load as autoload
urls = (
    "/hello", "hello"
)


class hello(object):
    def GET(self):
        return "hello"


if __name__ == "__main__":
    app = web.application(urls, globals())
    autoload.register(app)
    app.run()
