#!/bin/bash

PROG=../obj/bubble_ana
DIFF="diff -bBy"

EXPECT=vdetail.expect.dat
ACTUAL=vdetail.dat
DIFFOUT=diff.dat

$PROG  input.dat

if ! $DIFF $EXPECT $ACTUAL > $DIFFOUT
then
       printf "FAIL: Output Incorrect\n"
       printf "OUTPUT: EXPECT   vs   ACTUAL\n"
       cat $DIFFOUT
else
       printf "TEST RESULTS OK\n"
fi
