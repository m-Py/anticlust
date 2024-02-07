#!/bin/sh

R -e 'rmarkdown::render("README.Rmd")'
mv README.md ../README.md
cp man/figures/*.png ../man/figures
