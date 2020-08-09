#!/bin/sh
# make sure you have python and anaconda loaded

python MACIE_iter200.py -m predict -a ./output/exList_ready.txt -c ./parameters/Testing_Noncoding_class_prefix_f23.txt -o ./output/ -f 23 -p ./parameters/ -n exlist
