#!/bin/sh

R -e 'rmarkdown::render("README.Rmd")'
cp README.md ../README.md
cp man/figures/*.png ../man/figures
