# Theoretical Background

We search for the wave function of the ground state  in the form

![equation](https://latex.codecogs.com/gif.latex?%5CPsi%28%5Cvec%7B%5Ctheta%7D%29%20%3D%20U%28%5Cvec%7B%5Ctheta%7D%29%5CPsi_0)

where ![equation](https://latex.codecogs.com/gif.latex?%5CPsi_0) stands for the initial guess of the wave function and

![equation](https://latex.codecogs.com/gif.latex?U%28%5Cvec%7B%5Ctheta%7D%29%20%3D%20u_N%20%5Cdots%20u_1)

![equation](https://latex.codecogs.com/gif.latex?u_i%20%5Cequiv%20u%28%5Ctheta_i%29%20%3D%20e%5E%7B-i%5Ctheta_i%5Csigma/2%7D)

Equation on variation of the parameters

![equation](https://latex.codecogs.com/gif.latex?%5Csum_%7Bab%7D%20A_%7Bab%7D%20%5Cdot%7B%5Ctheta%7D_b%20%3D%20-C_a)

with

![equation](https://latex.codecogs.com/gif.latex?C_a%20%3D%20%5Cleft%5Clangle%20%5Cfrac%7B%5Cpartial%5Cpsi%7D%7B%5Cpartial%5Ctheta_a%7D%5Cleft%5Cvert%20H%5Cright%5Cvert%20%5Cpsi%20%5Cright%5Crangle%20&plus;%20%5Cleft%5Clangle%20%5Cpsi%20%5Cleft%5Cvert%20H%5Cright%5Cvert%20%5Cfrac%7B%5Cpartial%5Cpsi%7D%7B%5Cpartial%5Ctheta_a%7D%20%5Cright%5Crangle)

and

![equation](https://latex.codecogs.com/gif.latex?A_%7Bab%7D%20%3D%20%5Cleft%5Clangle%20%5Cfrac%7B%5Cpartial%5Cpsi%7D%7B%5Cpartial%5Ctheta_a%7D%20%5CBig%5Cvert%20%5Cfrac%7B%5Cpartial%5Cpsi%7D%7B%5Cpartial%5Ctheta_b%7D%20%5Cright%5Crangle%20&plus;%20%5Cleft%5Clangle%20%5Cfrac%7B%5Cpartial%5Cpsi%7D%7B%5Cpartial%5Ctheta_b%7D%20%5CBig%5Cvert%20%5Cfrac%7B%5Cpartial%5Cpsi%7D%7B%5Cpartial%5Ctheta_a%7D%20%5Cright%5Crangle)

Utilizing that the derivative with respect to the parameter expresses as

![equation](https://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Cpartial%20U%7D%7B%5Cpartial%20%5Ctheta_a%7D%20%3D%20-%5Cfrac%7Bi%7D%7B2%7Du_n%20%5Cdots%20u_a%20%5Csigma%5E%7B%28i_a%29%7D%20u_%7Ba-1%7D%20%5Cdots%20u_1%20%3D%20-%5Cfrac%7Bi%7D%7B2%7D%20V_a)

and decomposing the Hamiltonian into a sum of Pauli strings (PS)

![equation](https://latex.codecogs.com/gif.latex?H%20%3D%20%5Csum_h%20%5Calpha_h%20%5CSigma_h)

one can evaluate all elements of the equation. 
First, let us write PS as follows

![equation](https://latex.codecogs.com/gif.latex?%5CSigma_h%20%3D%20U_h%20D_h%20U%5E%5Cdagger_h)

where ![equation](https://latex.codecogs.com/gif.latex?D_h) is the diagonal matrix. Then the vector can be found as a result of the measurement of the circuit

![alt text](pictures/c_vec_circ.png)

as follows

![equation](https://latex.codecogs.com/gif.latex?C_a%20%3D%20%5Csum_x%20%5Cleft%5Clangle%20x%20%5Cleft%5Cvert%20D_h%20%5Cright%5Cvert%20x%20%5Cright%5Crangle%20%5Cleft%28%20P_%7B0x%7D%20-%20P_%7B1x%7D%20%5Cright%20%29)

# Notes

To launch the code type:
```shell
python prog.py z115_1-int.txt z115_2-int.txt 8
```
Here files z115_1-int.txt and z115_2-int.txt stores one- and two-body integrals in the specific format and number 8 stands for the number of spin-orbitals to be used.

The curent launching format will be changed soon and all necessary data: number of electrons, names of the files with one- and two-body integrals, number of orbitals, type of the anzatz, and etc. will be stored in a separate file.

# References 

S. McArdle, T. Jones, S. Endo, Y. Li, S.C. Benjamin, and Xiao Yuan,
``Variational ansatz-based quantum simulation of imaginary time evolution'',
[npj Quantum Inf. 5, 75 (2019)](https://www.nature.com/articles/s41534-019-0187-2)
