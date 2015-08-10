#!/bin/sh

HEADERS=`find include/ | grep \.h$ | grep -v logging`
SRCS=`find src/ | grep \.cxx$`

cppcheck --std=c++11 --enable=all $HEADERS $SRCS
