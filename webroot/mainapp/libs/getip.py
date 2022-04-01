# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

import web


def get_ip():
    ip = web.ctx.env.get('x-forwarded-for')
    if not ip or ip == "unknown":
        ip = web.ctx.env.get('x-forwarded-for')
    if not ip or ip == "unknonw":
        ip = web.ctx.env.get("Proxy-Client-IP")
    if not ip or ip == "WL-Proxy-Client-IP":
        ip = web.ctx.ip
    return ip
