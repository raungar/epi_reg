#!/bin/python

import sys
import urllib.request
import os

url=sys.argv[1]
file_name=sys.argv[2]

#if the file has not already been downloaded, please download
#done in this multi-script manner to fulfill urllib requirement
if not os.path.isfile(file_name):
	urllib.request.urlretrieve(url,file_name)
