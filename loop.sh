#!/bin/bash

set -e

while [ 1 ]
do
guid=`openssl rand -hex 32`
stdbuf -oL -eL ./traveling-santa 2>&1 > log_${guid}
done
