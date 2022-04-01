#!/bin/bash

convert $1".svg" -resize 600x600 $1".png"
convert $2".svg" -resize 600x600 $2".png"
convert $3".svg" -resize 600x600 $3".png"
convert $4".svg" -resize 600x600 $4".png"
