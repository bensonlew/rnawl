# -*- coding:utf-8 -*-
__author__ = 'guanqing.zou'
import re
import requests
from HTMLParser import HTMLParser
import os
import datetime
import random


proxy = [
    "114.239.149.127:808",
    "1.197.203.61:9999",
    "114.239.150.79:808",
    "1.198.72.85:9999",
    #"117.57.91.98:22413",
    "222.184.7.206:40908"
]



#登录源码：
def login(login_url,my_data, except_proxy=[]):
    #请求头
    my_headers = {
        'User-Agent' : 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_9_5) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/48.0.2564.116 Safari/537.36',
        'Accept' : 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
        'Accept-Encoding' : 'gzip',
        'Accept-Language' : 'zh-CN,zh;q=0.8,en;q=0.6,zh-TW;q=0.4',
        'Connection': 'keep-alive'
    }


    try_num = 1
    while try_num < 4:
        if except_proxy:
            new_proxy = list(set(proxy) - set(except_proxy))
        else:
            new_proxy = proxy

        try:
            if new_proxy:
                proxy_ip = random.choice(new_proxy)
                r = requests.post(login_url, headers = my_headers, data = my_data, proxies={'https': proxy_ip})
            else:
                r = requests.post(login_url, headers = my_headers, data = my_data)
            break
        except Exception as e:
            print(e)
            print("try num %s failed"% try_num)
            try_num +=1
            if new_proxy:
                print(proxy_ip)
                except_proxy.append(proxy_ip)
    if try_num == 4:
        raise  Exception("尝试3次连接仍失败")
    if r.status_code == 200:
        return  r.content
    else:
        raise Exception("return code is %s" % r.status_code)


class MyHTMLParser(HTMLParser):
    def __init__(self):
        HTMLParser.__init__(self)
        self.hit_path = ''
    def handle_starttag(self, tag, attrs):
        if tag == "a":
            if len(attrs) == 0: 
                pass
            else:
                for (variable, value) in attrs:
                    if variable == "href":
                        if 'gm_key_64.gz' in value:
                            self.hit_path = value
    

def check_data(log_dir,update_day=30):
    update = False
    gm_key_log = os.path.join(log_dir, '.gm_key_log')
    current_date = datetime.datetime.now()
    current_date_str = current_date.strftime("%Y-%m-%d %H:%M:%S")
    if not os.path.exists(gm_key_log):
        print('Not Found .gm_key_log')
        update = True
    else:
        with open(gm_key_log) as fr:
            data_log = fr.readlines()
            last_log = data_log[-1].strip()
            if last_log:
                try:
                    pre_data = datetime.datetime.strptime(last_log, "%Y-%m-%d %H:%M:%S")

                    diff_day = (current_date - pre_data).days
                    if diff_day > update_day:
                        print('Over Day %s'%diff_day)
                        update = True
                except Exception as e:
                    print('Read log WRONG')
                    update = True
            else:
                update = True

    return update,current_date_str, gm_key_log


def update_gm_key_pip(result_dir='./', rely_biocluster=True):
    if rely_biocluster:
        from biocluster.config import Config
        config = Config()
        home_path = os.path.dirname(config.SOFTWARE_DIR)
    else:
        home_path = os.environ['HOME']
    update, current_data, log_file = check_data(home_path, 30)
    if update:
        my_data = {
            'name': 'major',
            'institution': 'major',
            'address':'',
            'city': '',
            'state': '',
            'country': 'china',
            'email': 'zouguanqing@majorbio.com',
            'submit': 'I agree to the terms of this license agreement',
            'os': 'linux64',
            'program': 'gms'
        }

        url= 'http://topaz.gatech.edu/GeneMark/license_download.cgi'
        content = login(url, my_data)
        hp = MyHTMLParser()
        hp.feed(content)
        hp.close()
        print(hp.hit_path)
        key_obj=requests.get(hp.hit_path)
        result_file = os.path.join(result_dir, 'gm_key_64.gz')
        decompress_result = os.path.join(result_dir, 'gm_key_64')
        if os.path.exists(result_file):
            os.remove(result_file)
        with open(result_file,'w') as fw:
            fw.write(key_obj.content)
        os.system('zcat {} > {}'.format(result_file, decompress_result))
        os.system('cp {} {}/.gm_key'.format(decompress_result, home_path))
        print('Update Done')
        with open(log_file, 'a') as flog:
            flog.write('\n%s'%current_data)

        print('Update %s'%log_file)
        restr = 'Update Done'

    else:
        restr = 'NOT Update'
        print('NOT Update')

    return restr



if __name__ == '__main__':
    update_gm_key_pip()








