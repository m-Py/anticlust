
#!/bin/sh

R -e "source('build-package.R')"

# copy vignette to my website
cp ./inst/doc/my-vignette.html ../m-Py.github.io/anticlust/stimulus-selection.html

cd ../m-Py.github.io/anticlust

git add stimulus-selection.html
#git commit -m "added latest version of vignette"

cd -
