#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import string
import shutil

toDir = 'yale2/';
fromDir = 'yaleB18/';

files=os.listdir(fromDir);
files.sort();
counter=1;
for j in range(len(files)):
	src = fromDir + '/' + files[j];
	dst = toDir + str(counter) + '.pgm'
	shutil.copyfile(src, dst)
	counter = counter + 1;


