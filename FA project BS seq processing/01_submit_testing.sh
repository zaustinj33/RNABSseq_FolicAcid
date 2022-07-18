#!/bin/bash

# Usage: ./01_submit.sh <File with list of data> <Absolute path to directory>

cat $1 | while read LINE
do
	sh ./07_CallCutoff_merged.sh $LINE $2
done
