# Ewald_Summation

### Members: [Zhe Chen](zc1291@cims.nyu.edu); [Guanchun Li](guanchun.li@nyu.edu)

HPC Final Project, 2019 spring, Courant, NYU

---

### Fast parallel cuda code for Ewald Summation for Skokes potential

I, with Guanchun Li, will together take the Ewald Summation problem for our final prokect for HPC class. Here's detail of our plan.

Green function of Stokes flow is shown in Oseen-Burgers tensor, or Stokeslet. Thus, convolution of green function is the key to compute velocity field. In free space, the main difficulty is the slow decay, 1/r, of the kernel, which requires big N and O(N^2) complexity. Fortunately, Ewald Summation method does great jobs in this problem. Generally, it split the convolution into real space, which can be adjusted by Ewald parameter $\xi$, and k-space, which can be computed fast by FFT based method. By  this mean, Ewald method change convolution of this problem from O(N^2) to O(N logN) with spectral accuracy.

The fact that Ewald method is FFT based makes it possible for us to do parallel computing. We plan to write cuda gpu-parallel code to improve Ewald method. First, we want to implemented parallel FFT method by our own or using cuFFT library, which can be plugged into solution to k-space of  Ewald decomposition.  Also, there's an integral step in algorithm of Ewald Summation, which is described in <http://dx.doi.org/10.1016/j.jcp.2010.08.026>. Since it's periodic, we can simply use trapezoidal integral to get spectral accuracy. Moreover, the trapezoidal integral can also be paralleled since it's a summation of very long vector.

The goal of this project is to develop a CUDA-based gpu-parallel libary to implement Ewald summation fast. Hopefully, we could use my CFD final project, which is about particles in a stokes flow above a wall, as an exemple to discuss  performance of this method.


### GPU library for NuFFT
https://github.com/andyschwarzl/gpuNUFFT/tree/master/CUDA

### Usage:
Environment on cims machine:

1) cuda-8.0

2) gcc-4.9              

3) mpi/openmpi-x86_64

```
cd src/
make
./main num_threads N num_p P rp
```
