# -*- coding: utf-8 -*-
# __author__ = 'fwy'
import cairosvg
import os
import argparse

def svg2pdf(svg_path,pdf_path):
    cairosvg.svg2pdf(url=svg_path, write_to=pdf_path)

def svg2png(svg_path,png_path):
    cairosvg.svg2png(url=svg_path, write_to=png_path)