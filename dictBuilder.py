#!/bin/python3

import codecs
import sys
import dbm

if len(sys.argv) >= 2:
    file = open(sys.argv[1], 'rb')
    weather = dbm.open('weather', 'c')
    for line in file:
        content = codecs.decode(line, encoding='utf-8')
        key, value =  content.strip('\r\n').split(',')
        weather[key] = value

    file.close()
    weather.close()
else:
    sys.exit(1)
