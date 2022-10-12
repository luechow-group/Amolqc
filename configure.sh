#!/bin/bash

if [[ -e "$1" ]]; then
	target="$1"
elif [[ -e "cmake/configs/$1" ]]; then
    target="cmake/configs/$1"
elif [[ -e "cmake/configs/$1.cmake" ]]; then
    target="cmake/configs/$1.cmake"
else
	echo "Config file '$1' not found"
	exit 1
fi

cp "$target" cmake/config.cmake
echo "Copied $target to cmake/config.cmake."
