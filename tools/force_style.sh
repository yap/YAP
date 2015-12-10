#!/bin/sh

astyle --options=tools/astylerc --formatted \
    include/*.h \
    src/*.cxx \
    programs/*.cxx \
    examples/*/*.cxx \
    examples/*/*.h \
    test/*.cxx