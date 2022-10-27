from distutils.core import setup, Extension
import os
dir = os.getcwd()
module1 = Extension('VDmod',
                    sources = ['VDmodmodule.c'],
                    libraries = ['VDIntegral'],
                    library_dirs = [dir],
                    include_dirs = [dir]
                    )

setup (name = 'VDmod',
       version = '1.0',
       description = 'Test to solve integral VD',
       ext_modules = [module1])

#python setup.py build
#or
#python setup.py install