#!/usr/bin/python
# coding: utf-8
# __author__ = 'Zhao Binbin'
# 自动发邮件模块

import os
import re
import smtplib
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText


class SendEmail(object):
    """
    自动发邮件脚本
    """

    def __init__(self, sender, smtpserver, password, username, receiver, subject, port=25):
        self.sender = sender
        self.smtpserver = smtpserver
        self.password = password
        self.username = username
        self.receiver = receiver
        self.subject = subject
        self.msg = MIMEMultipart('mixed')
        self.msg['Subject'] = subject
        self.msg['From'] = sender
        self.msg['To'] = ";".join(receiver.strip().split(","))
        self.port = port

    def send_msg(self, content):
        """
        content 为要书写的内容
        :param content:
        :return:
        """
        content = content
        content_plain = MIMEText(content, 'plain', 'utf-8')
        self.msg.attach(content_plain)

    def send_image(self, file_path, file_name=None):
        """
        file_path输入要上传图片的路径
        :param file_path:
        :param file_name:
        :return:
        """
        if file_name is not None:
            file_name = file_name
        else:
            file_name = os.path.basename(file_path)
        sendimagefile = open(file_path, 'rb').read()
        image = MIMEImage(sendimagefile)
        image.add_header('Content-ID', '<image1>')
        image.add_header('Content-Disposition', 'attachment', filename=file_name)
        self.msg.attach(image)

    def send_html(self, link, file_name=None):
        """
        link 为html格式，举例如下：
        html = <a href="http://www.baidu.com">link</a>
        :param link:
        :param file_name:
        :return:
        """
        html = link
        if file_name is not None:
            file_name = file_name
        else:
            result = "".join(re.findall(".*http(.*)com.*", html))
            file_name = "http" + result + "com"
        text_html = MIMEText(html, 'html', 'utf-8')
        text_html.add_header('Content-Disposition', 'attachment', filename=file_name)
        self.msg.attach(text_html)

    def attachment(self, file_path, file_name=None):
        """
        file_path为附件地址
        :param file_path:
        :param file_name:
        :return:
        """
        sendfile = open(file_path, 'rb').read()
        text_att = MIMEText(sendfile, 'base64', 'utf-8')
        text_att["Content-Type"] = 'application/octet-stream'
        if file_name is not None:
            file_name = file_name
        else:
            file_name = os.path.basename(file_path)
        text_att.add_header('Content-Disposition', 'attachment', filename=file_name)
        self.msg.attach(text_att)

    def send_email(self):
        server = smtplib.SMTP_SSL(self.smtpserver, self.port)  # 发件人邮箱中的SMTP服务器，端口是25
        server.login(self.username, self.password)  # 括号中对应的是发件人邮箱账号、邮箱密码
        server.sendmail(self.sender, self.receiver, self.msg.as_string())  # 括号中对应的是发件人邮箱账号、
        # 收件人邮箱账号、发送邮件
        server.quit()  # 关闭连接


if __name__ == "__main__":
    a = SendEmail("1274095891@qq.com", "smtp.qq.com", "ocfnjnhnsbqubaej", "1274095891@qq.com",
                  "1274095891@qq.com,bbzhao@ips.ac.cn,nongwuyunfeng@163.com",
                  "来自20181008下午邮件的自动发送测试", 465)
    a.send_msg("我是中国人")
    a.send_image("C:\\Users\\binbin.zhao\\Desktop\\BAD810EA-C119-432c-B06A-75F5D621C2F4.png")
    a.attachment("C:\\Users\\binbin.zhao\\Desktop\\content of mail format.txt")
    a.send_html("<a href=\"http://www.baidu.com\">link</a>")
    a.send_email()
