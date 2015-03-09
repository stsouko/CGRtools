#!/bin/bash

if [ ! -d "./pyinstaller" ] ; then
  echo "download https://github.com/pyinstaller/pyinstaller into ./pyinstaller dir"
  exit 1
fi

if [ ! -d "./pyinstaller" ] ; then
  echo "download http://upx.sourceforge.net into ./upx-lin dir"
  exit 1
fi

python ./pyinstaller/pyinstaller.py --upx-dir=./upx-lin -y -F condenser.py

