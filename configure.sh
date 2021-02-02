#!/bin/bash

if [[ -e "$1" ]]; then
	target="$1"
elif [[ -e "make/configs/$1" ]]; then
	target="make/configs/$1"
elif [[ -e "make/configs/$1.mk" ]]; then
    target="make/configs/$1.mk"
    if [[ -e "cmake/configs/$1.cmake" ]]; then
        target2="cmake/configs/$1.cmake"
    fi
elif [[ -e "cmake/configs/$1" ]]; then
    target="cmake/configs/$1"
elif [[ -e "cmake/configs/$1.cmake" ]]; then
    target="cmake/configs/$1.cmake"
else
	echo "Config file '$1' not found"
	exit 1
fi

if [[ $target =~ "cmake" ]]; then
    cp "$target" cmake/config.cmake
    echo "Copied $target to cmake/config.cmake."
else
    cp "$target" make/config.mk
    echo "Copied $target to make/config.mk."
fi

if [[ ! -z "$target2" ]]; then
    cp "$target2" cmake/config.cmake
    echo "Copied $target2 to cmake/config.cmake."
fi
