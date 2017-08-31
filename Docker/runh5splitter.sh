#!/bin/bash

Rscript --vanilla --default-packages=methods,stats,utils /usr/local/bin/runh5split4batman.R -i $1 -t $2
