#coding:utf-8
import requests
import os
from functools import wraps
from tqdm import tqdm

def try_for_good(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except:
            return wrapper(*args, **kwargs)
    return wrapper

@try_for_good
def download(filepath=None, save_name=None):
    email= "zhangyujuan@webmail.hzau.edu.cn"
    password = "Majorbio@123"
    url= "https://cancer.sanger.ac.uk/cosmic/file_download/"
    r = requests.get(url+filepath, auth=(email, password))
    download_url = r.json()["url"]
    ## 支持断点续传
    response = requests.get(download_url, stream=True)
    file_size = int(response.headers['content-length'])
    if os.path.exists(save_name):
        print("Resume from break-point: {}".format(download_url))
        first_byte = os.path.getsize(save_name)
    else:
        first_byte = 0
        print("Downloading: {}".format(download_url))
    header = {"Range": "bytes={}-{}".format(first_byte, file_size)}
    pbar = tqdm(total=file_size, initial=first_byte, unit='B', unit_scale=True, desc=save_name)
    link = requests.get(download_url, headers=header, stream=True)
    with open(save_name, "ab") as f:
        for chunk in link.iter_content(chunk_size=10240):
            if chunk:
                f.write(chunk)
                pbar.update(10240)
    pbar.close()
    cmd = "gunzip {}".format(save_name)
    print(cmd)
    code = os.system(cmd)
    if code == 0:
        print("Successed downloading: {}".format(save_name))
        return True
    else:
        raise Exception("Failed downloading: {}".format(save_name))
    

def main():
    # get a list of all files in COSMIC v91
    filelist = requests.get("https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v91")

    # extract the download URLs from the JSON response
    for filepath in filelist.json():
    # extract the filename from the download file path
        filename = filepath.rsplit("/", 1)[1]
        if os.path.exists(filename.split(".gz")[0]):
            print("skip download file {}".format(filename))
            continue
    # skip the Oracle database dump; it's around 50Gb
        if "ORACLE_EXPORT" in filename:
            print("skipping oracle export file")
            continue
    # get the URL for downloading the file
        link = download(filepath=filepath, save_name=filename)
        
    
if __name__ == "__main__":
    main()
    
