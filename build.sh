#!/bin/bash
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ ! -d "./pyinstaller" ] ; then
    echo "downloading pyinstaller into ./pyinstaller"
    git clone git@github.com:pyinstaller/pyinstaller.git
fi

if [ ! -d "./upx-lin" ] ; then
    echo "downloading upx into ./upx-lin"
    wget http://upx.sourceforge.net/download/upx-3.91-amd64_linux.tar.bz2 -O - | tar -xj
    mv upx-3.91-amd64_linux upx-lin
fi

nuitka --standalone --portable --python-version=3.4 --recurse-all --recurse-stdlib main.py

python3 ./pyinstaller/pyinstaller.py  --additional-hooks-dir=./hooks --upx-dir=./upx-lin --name=cgrtools --key=jyvt0n3 -y -F main.py

