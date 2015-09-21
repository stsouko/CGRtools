#!/bin/bash
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ ! -d "./pyinstaller" ] ; then
    echo "downloading pyinstaller into ./pyinstaller"
    git clone git@github.com:pyinstaller/pyinstaller.git
fi

python3 ./pyinstaller/pyinstaller.py  --additional-hooks-dir=./hooks --name=cgrtools --key=$RANDOM$RANDOM$RANDOM -y -F main.py
