#-*- coding:utf-8 -*-
from datetime import datetime
import requests
import bs4
import pandas as pd
import sys
import time
reload(sys)
sys.setdefaultencoding("utf-8")

def main(name_type, plant_list, output_path):

    df_list = list()
    for i in plant_list.split(','):
        time.sleep(1)
        with requests.Session() as session:
            headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/52.0.2743.116 Safari/537.36"}
            headers_ = {
              "Host":"pe.ibcas.ac.cn",
              "User-Agent":"Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:88.0) Gecko/20100101 Firefox/88.0",
              "Accept":"text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
              "Accept-Language":"zh-CN,zh;q=0.8,zh-TW;q=0.7,zh-HK;q=0.5,en-US;q=0.3,en;q=0.2",
              "Accept-Encoding":"gzip, deflate",
              "Content-Type":"application/x-www-form-urlencoded",
              "Content-Length":"374",
              "Origin":"http://pe.ibcas.ac.cn",
              "Connection":"keep-alive",
              "Referer":"http://pe.ibcas.ac.cn/herbs/plantnames.aspx",
              "Upgrade-Insecure-Requests":"1"
            }
            r = session.get('http://pe.ibcas.ac.cn/herbs/plantnames.aspx', headers=headers)
            soup = bs4.BeautifulSoup(r.content, "html5lib")
            viewState = soup.select_one('#__VIEWSTATE')["value"]
            viewStateGenerator = soup.select_one('#__VIEWSTATEGENERATOR')["value"]
            eventValidation = soup.select_one('#__EVENTVALIDATION')["value"]
            if name_type == 'latin':
                paramsPost = {'__VIEWSTATE': viewState,
                        '__VIEWSTATEGENERATOR': viewStateGenerator,
                        '__EVENTVALIDATION': eventValidation,
                        "TextBox1":i,
                        "TextBox2":"",
                        "TextBox3":"",
                        "ImageButton1.x":"43",
                        "ImageButton1.y":"9"
                        }
            elif name_type == 'chinese':
                paramsPost = {'__VIEWSTATE': viewState,
                        '__VIEWSTATEGENERATOR': viewStateGenerator,
                        '__EVENTVALIDATION': eventValidation,
                        "TextBox1":"",
                        "TextBox2":i,
                        "TextBox3":"",
                        "ImageButton1.x":"43",
                        "ImageButton1.y":"9"
                        }
            elif name_type == 'english':
                paramsPost = {'__VIEWSTATE': viewState,
                        '__VIEWSTATEGENERATOR': viewStateGenerator,
                        '__EVENTVALIDATION': eventValidation,
                        "TextBox1":"",
                        "TextBox2":"",
                        "TextBox3":i,
                        "ImageButton1.x":"43",
                        "ImageButton1.y":"9"
                        }
            response = session.post("http://pe.ibcas.ac.cn/herbs/plantnames.aspx", data=paramsPost, headers=headers_)
            if response.status_code == 200:
                response_ = session.get("http://pe.ibcas.ac.cn/herbs/plantnameres.aspx")
                soup = bs4.BeautifulSoup(response_.content, "html5lib")
                # print(soup)
                table = soup.find_all(id='DataGrid1')[0]
                df = pd.read_html(str(table),encoding="utf-8")[0]
                df.drop([0],inplace=True)
                df_list.append(df)
                # df.to_csv('plant.txt', header=None, index=None, sep='\t')
    all_df = reduce(lambda x, y: pd.concat([x, y], axis=0), df_list)
    all_df.fillna("")
    all_df.columns = ['Latin name', 'Chinese name', 'English name']
    all_df.to_csv(output_path, header=True, index=None, sep='\t')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", type=str, help="name type")
    parser.add_argument("-p", type=str, help="plant list")
    parser.add_argument("-o", type=str, help="output path")
    args = parser.parse_args()
    main(args.n, args.p, args.o)
