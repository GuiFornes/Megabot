ARDUINO_DIR  = /usr/share/arduino
USER_LIB_PATH=.
ARDUINO_LIBS = LiquidCrystal-I2C Wire
BOARD_TAG    = uno
ARDUINO_PORT?= /dev/ttyUSB0

CONTROLER_ID?=0
CFLAGS=-DID=$(CONTROLER_ID)
CXXFLAGS=-DID=$(CONTROLER_ID)
AVRDUDE_OPTS=-v

#avoid compiling simulator:
LOCAL_CPP_SRCS=

include /usr/share/arduino/Arduino.mk
