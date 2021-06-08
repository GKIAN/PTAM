# PTAM

Peak Trough Averaging Method (峰谷平均法).

## License

[The MIT License](http://tchel.mit-license.org/)

## Author

Tche LIU, seistche@gmail.com, USTC

## Description

[mod_ptam](mod_ptam.F90) is a `fortran` module used to implement the peak trough averaging method to an external function.

To illustrate how to use the fortran module, I take the equation (6b) in (张海明等, 2001) as an example:

$ \int_0^{+ \infty} e^{ - \alpha k} J_1 (\beta k) dk = \frac{1}{\beta} \left[ 1 - \frac{\alpha}{\sqrt{\alpha^2 + \beta^2}} \right] $

In the source file [main.F90](main.F90), I defined a function `exfunc` followed the example formula. After compiling the source files by `gfortran mod_ptam.F90 main.F90 -o a.out`, we can run the example by `./a.out`, and the output will be similar with the following:

```shell
The PTAM result is: 0.333209425
The analytical result is: 0.333222240
The relative error is: -0.384578125E-4
```

## Reference

- 张海明, 陈晓非, 张似洪, (2001). 峰谷平均法及其在计算浅源合成地震图中的应用. _地球物理学报_: 44(6), 805-813
- Zhang, H. M. ,  Chen, X. F. , &  Chang, S., (2003). An efficient  numerical method for computing synthetic seismograms for a layered  half-space with sources and receivers at close or same depths. _Birkhäuser Basel_.