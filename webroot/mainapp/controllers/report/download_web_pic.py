# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
# from gevent import monkey; monkey.patch_all()
from mainapp.libs.signature import check_sig
import web
import json
import tempfile
import HTMLParser
from subprocess import Popen
# from gevent.subprocess import Popen
from biocluster.api.file.remote import RemoteFileManager
import random
import datetime
from biocluster.config import Config
import os

if "LD_LIBRARY_PATH" not in os.environ.keys():
    os.environ['LD_LIBRARY_PATH'] = ""
os.environ['LD_LIBRARY_PATH'] = Config().SOFTWARE_DIR + '/program/Python35/lib:' + os.environ['LD_LIBRARY_PATH']


class DownloadWebPic(object):
    """
    网页的svg转化成pdf或者png等图片格式
    POST参数有(均为必须参数) "file_type", "file_name", "svg_data", "scale"
    分别为 文件类型(必须为png或者pdf), 文件名称(暂时没有用处，但是会放在header中),svg数据,放大大小(在png模式下生效)
    """

    @check_sig
    def POST(self, my_type=None):
        data = web.input()
        print('Convert Type: {}'.format(my_type))
        for i in ["file_type", "file_name", "svg_data", "scale"]:
            if not hasattr(data, i):
                msg = {"success": False, "info": "缺少参数: {}".format(i)}
                return json.dumps(msg)
            if getattr(data, i).strip() == "":
                msg = {"success": False, "info": "参数%s不能为空" % i}
                return json.dumps(msg)
        if data.file_type not in ['png', 'pdf']:
            msg = {"success": False, "info": "结果文件类型必须为png或者pdf:{}".format(data.file_type)}
            return json.dumps(msg)
        if not data.scale.isalnum():
            msg = {"success": False, "info": "scale参数必须是数值:{}".format(data.scale)}
            return json.dumps(msg)
        file_pic = self._svg_convert()
        if file_pic:
            # web.header('Content-Type', 'application/octet-stream')
            # web.header('Transfer-Encoding', 'chunked')
            # web.header('Access-Control-Allow-Origin', '*')
            # web.header('Content-disposition', 'attachment; filename={}'.format(data.file_name))
            # return open(file_pic, 'rb').read()
            msg = {'success': True, 'info': file_pic}
        else:
            msg = {"success": False, "info": "生成图片文件出错".format(data.scale)}
        print(json.dumps(msg))
        return json.dumps(msg)

    def _svg_convert(self):
        """
        转换svg为图片，使用python3的包cairosvg命令行形式生成
        """
        temp_dir = tempfile.mkdtemp()
        temp_svg = temp_dir + '/temp.svg'
        with open(temp_svg, 'wb') as w:
            parser = HTMLParser.HTMLParser()
            svg_data = parser.unescape(web.input().svg_data)  # html转义
            w.write(svg_data)
        file_name = random_file_name() + '.' + web.input().file_type
        temp_pic = temp_dir + '/' + file_name
        cmd = Config().SOFTWARE_DIR + '/program/Python35/bin/'
        cmd += 'cairosvg {} -f {} -o {} -s {}'.format(temp_svg, web.input().file_type, temp_pic, int(web.input().scale))
        pro = Popen(cmd, shell=True)
        pro.wait()
        if pro.returncode == 0:
            print("TEMP PIC: {}".format(temp_pic))
            target_dir = self.upload_pic(temp_pic)
            print "return target_dir : %s" % target_dir
            if target_dir:
                return target_dir + file_name
            else:
                return
        else:
            return

    def upload_pic(self, pic):
        """
        上传结果文件
        """
        target = self.create_remote_target()
        remote = RemoteFileManager(target)
        try:
            remote._type = "local"
            remote.upload(pic)
        except Exception as e:
            print('UPLOAD ERROR: {}'.format(e))
            return
        else:
            return target

    def create_remote_target(self):
        """
        构建远程目录结构/位置
        """
        data = web.input()
        target_dir = 'tsanger:'
        target_dir = 's3://'
        if hasattr(data, 'client') and data.client == 'client01':
            target_dir = 'sanger:'
            target_dir = 's3nb://'
        target_dir += 'rerewrweset/report_img/download_img/'
        # file_name = random_file_name()
        # target = target_dir
        print("REMOTE DIR: {}".format(target_dir))
        return target_dir


def random_file_name():
    sed = 'qwertyuiopasdfghjklzxcvbnm1234567890'
    random_value = ''.join(random.sample(sed, 4))
    date = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    return date + random_value



if __name__ == '__main__':
    pass
