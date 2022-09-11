#!/bin/bash

base_repo=$1
name=$2
ranking_file=$3
compress_file=$4

if test -f "$base_repo$name$ranking_file"; then
	touch $base_repo$name$ranking_file
	
	if test -f "$base_repo$name"/*"$compress_file"; then
		touch $base_repo$name/*$compress_file
	fi

else
	if test -f "$base_repo$name"/*"$compress_file"; then
		echo "$base_repo$name$ranking_file does not exist\t removing $base_repo$name/*$compress_file\n"
		rm $base_repo$name/*$compress_file
	fi

fi
