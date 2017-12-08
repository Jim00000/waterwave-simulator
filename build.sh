nvcc -Xcompiler -fPIC -O2 -c wave_equation_cuda_api.cu
python setup.py build_ext --inplace