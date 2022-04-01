# -*- coding: utf-8 -*-

import os
import json
import logging
import re
import sys
import regex

from biocluster.api.file.remote import RemoteFileManager
from biocluster.api.file.lib.transfer import TransferManager
import glob
from biocluster.api.file.lib.transfer import MultiFileTransfer


def download_from_s3(from_file, to_path, cover=True):
    if not to_path.startswith("/"):
        to_path = os.path.join(os.getcwd(), to_path)
        print to_path
    if os.path.exists(to_path):
        os.remove(to_path)
    transfer = MultiFileTransfer()
    transfer.add_download(from_file, to_path, base_path=from_file)
    transfer.perform()

if __name__ == '__main__':
    from_path = sys.argv[1]
    to_path = sys.argv[2]
    download_from_s3(from_path, to_path)
    

