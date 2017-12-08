"""
Copyright (C) 2017 the team of Jim00000, ActKz and pityYo

Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
and associated documentation files (the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial 
portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT 
LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

setup(
    ext_modules = cythonize(Extension(
        "wave_equation",                                                # the extension name
        sources=["wave_equation.pyx", "wave_equation_api.cpp"],         # the Cython source and additional C source files
        language="c++",                                                 # the language
        libraries=["m"],                                                # linking libraries
        extra_compile_args=["-std=c++11", "-O2", "-Wall", "-Wextra", 
        "-fopenmp"],                                                    # compiler flags
        extra_link_args=["-fopenmp", "wave_equation_cuda_api.o",
        "-lcuda", "-lcudart", "-lboost_system", "-lboost_thread"],      # link arguments
        include_dirs=[numpy.get_include()],                             # included directories
)))