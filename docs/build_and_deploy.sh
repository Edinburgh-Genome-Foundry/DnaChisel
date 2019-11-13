#!/usr/bin/env bash

make html
cd _build/html
git add .
git commit -m "New docs"
git push origin gh-pages
cd ../..