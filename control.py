#!/usr/bin/python


import serial
import sys

serial.Serial(sys.argv[1], 115200, timeout=.1)



