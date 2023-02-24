from setuptools import Extension, setup

module = Extension("mykmeanssp", sources=['spkmeansmodule.c', 'matrix.c'])
setup(name='mykmeanssp',
     version='1.0',
     description='Python wrapper for the implementation of the spk, wam, ddg, gl and jacobi algorithms',
     ext_modules=[module])