#!/bin/bash
cd /usr/local/bin
Rscript --vanilla --default-packages=methods,stats,utils runh5split4batman.R -i $1 -t $2
