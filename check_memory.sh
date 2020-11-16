#!/bin/bash

R -e "devtools::install()"

R -d valgrind --vanilla -e "source('testing_valgrind.R')"
