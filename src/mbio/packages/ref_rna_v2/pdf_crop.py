# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
import sys
import os
from pdfCropMargins import crop

if __name__ == '__main__':
    pdf_in = sys.argv[1]
    pdf_out = sys.argv[2]

    crop(["-p", "1",  pdf_in, "-o", pdf_in + ".tmp.pdf"])
    crop(["-a", "-100", pdf_in + ".tmp.pdf", "-o", pdf_out])
    os.remove(pdf_in + ".tmp.pdf")
