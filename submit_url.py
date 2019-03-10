#!/bin/python

import sys
import urllib.request
import os

url=sys.argv[1]
file_name=sys.argv[2]


if not os.path.isfile(file_name):
	urllib.request.urlretrieve(url,file_name)
