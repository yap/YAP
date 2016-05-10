#!/bin/sh

astyle --options=tools/astylerc --formatted \
    include/*.h \
    src/*.cxx \
    examples/*/*.cxx \
    examples/*/*.h \
    test/*.cxx
