from distutils.core import setup

setup(
    name='condenser',
    version='1.0',
    packages=['condenserpkg'],
    url='https://github.com/stsouko/condenser',
    license='GPLv2',
    author='stsouko',
    author_email='stsouko@live.ru',
    description='simple CLI CGR generator',
    scripts=['condenser.py'], requires=['numpy']
)
