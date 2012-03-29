#!/bin/bash

echo "$@"
echo $@
echo $0
echo $1
echo $2


args=("$@")
fcID=${args[0]}
project_id=${args[1]}
bedfile=${args[2]}

