nvcc -O3 main.cu -o cuda_main -lcufft -arch=sm_61 -Xcompiler "/wd 4819"