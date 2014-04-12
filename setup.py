from distutils.core import setup
from condenserpkg.version import version

setup(
    name='condenser',
    version=version(),
    packages=['condenserpkg'],
    url='https://github.com/stsouko/condenser',
    license='GPLv2',
    author='stsouko',
    author_email='stsouko@live.ru',
    description='simple CLI CGR generator',
    scripts=['condenser.py'], requires=['numpy']
)
