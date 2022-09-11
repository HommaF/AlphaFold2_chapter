#!/bin/bash

repo=$1
name=$2
file=$3


if test -f "$repo$name$file"; then
	touch $repo$name$file

else
	echo "$repo$name$file does not exist"
fi
