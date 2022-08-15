#!/usr/bin/env bash
echo "self, scpt=${0}"; scpt=${0}

docker run -it --rm -v $PWD:$PWD -w $PWD rocker/tidyverse bash
