nvcc -Xcompiler -fPIC -O2 -c wave_equation_cuda_api.cu -D_FORCE_INLINES
python setup.py build_ext --inplace