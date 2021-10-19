# Theoretical Background

We search for the wave function of the ground state  in the form

![equation](https://latex.codecogs.com/gif.latex?%5CPsi%28%5Cvec%7B%5Ctheta%7D%29%20%3D%20U%28%5Cvec%7B%5Ctheta%7D%29%5CPsi_0)

where ![equation](https://latex.codecogs.com/gif.latex?%5CPsi_0) stands for the initial guess of the wave function and

![equation](https://latex.codecogs.com/gif.latex?U%28%5Cvec%7B%5Ctheta%7D%29%20%3D%20u_N%20%5Cdots%20u_1)

In the simplest case, ![equation](https://latex.codecogs.com/gif.latex?u_i) designates the one qubit rotation and expresses as

![equation](https://latex.codecogs.com/gif.latex?u_i%20%5Cequiv%20u%28%5Ctheta_i%29%20%3D%20e%5E%7B-i%5Ctheta_i%5Csigma/2%7D)

In order to find the parameters ![equation](https://latex.codecogs.com/gif.latex?%5Cvec%7B%5Ctheta%7D) one can use the imaginary time evolution (ITE) method.
Within this approach, the application of the ![equation](https://latex.codecogs.com/gif.latex?e%5E%7B-%5Ctau%20H%7D) to the wave function ![equation](https://latex.codecogs.com/gif.latex?%5CPsi%28%5Cvec%7B%5Ctheta%7D_%7B%5Crm%20old%7D%29) is equivalent to the following change of the parameters 

![equation](https://latex.codecogs.com/gif.latex?%5Cvec%7B%5Ctheta%7D_%7B%5Crm%20old%7D%20%5Crightarrow%20%5Cvec%7B%5Ctheta%7D_%7B%5Crm%20new%7D%20%3D%20%5Cvec%7B%5Ctheta%7D_%7B%5Crm%20old%7D%20&plus;%20%5Ctau%20%5Cdot%7B%5Cvec%7B%5Ctheta%7D%7D)

where ![equation](https://latex.codecogs.com/gif.latex?%5Cdot%7B%5Cvec%7B%5Ctheta%7D%7D) is found from the equation

![equation](https://latex.codecogs.com/gif.latex?%5Csum_%7Bb%7D%20A_%7Bab%7D%20%5Cdot%7B%5Ctheta%7D_b%20%3D%20-C_a)

Here

![equation](https://latex.codecogs.com/gif.latex?C_a%20%3D%20%5Cleft%5Clangle%20%5Cfrac%7B%5Cpartial%5Cpsi%7D%7B%5Cpartial%5Ctheta_a%7D%5Cleft%5Cvert%20H%5Cright%5Cvert%20%5Cpsi%20%5Cright%5Crangle%20&plus;%20%5Cleft%5Clangle%20%5Cpsi%20%5Cleft%5Cvert%20H%5Cright%5Cvert%20%5Cfrac%7B%5Cpartial%5Cpsi%7D%7B%5Cpartial%5Ctheta_a%7D%20%5Cright%5Crangle)

and

![equation](https://latex.codecogs.com/gif.latex?A_%7Bab%7D%20%3D%20%5Cleft%5Clangle%20%5Cfrac%7B%5Cpartial%5Cpsi%7D%7B%5Cpartial%5Ctheta_a%7D%20%5CBig%5Cvert%20%5Cfrac%7B%5Cpartial%5Cpsi%7D%7B%5Cpartial%5Ctheta_b%7D%20%5Cright%5Crangle%20&plus;%20%5Cleft%5Clangle%20%5Cfrac%7B%5Cpartial%5Cpsi%7D%7B%5Cpartial%5Ctheta_b%7D%20%5CBig%5Cvert%20%5Cfrac%7B%5Cpartial%5Cpsi%7D%7B%5Cpartial%5Ctheta_a%7D%20%5Cright%5Crangle)

The coefficients ![equation](https://latex.codecogs.com/gif.latex?A_%7Bab%7D) and ![equation](https://latex.codecogs.com/gif.latex?C_a) can be evaluated with the use of a quantum computer as follows.
First, let us decompose the Hamiltonian into a sum of Pauli strings (PS)

![equation](https://latex.codecogs.com/gif.latex?H%20%3D%20%5Csum_h%20%5Calpha_h%20%5CSigma_h)

As a next step, one needs to evaluate the derivative with respect to the parameter ![equation](https://latex.codecogs.com/gif.latex?%5Ctheta_a).
In the simplest case, when ![equation](https://latex.codecogs.com/gif.latex?u_i) stands for the one qubit rotations, one obtains

![equation](https://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Cpartial%20U%7D%7B%5Cpartial%20%5Ctheta_a%7D%20%3D%20-%5Cfrac%7Bi%7D%7B2%7Du_n%20%5Cdots%20u_a%20%5Csigma%5E%7B%28i_a%29%7D%20u_%7Ba-1%7D%20%5Cdots%20u_1%20%3D%20-%5Cfrac%7Bi%7D%7B2%7D%20V_a)

Utilizing this expression and the decomposition of the Hamiltonian into PS, one can show that

![equation](https://latex.codecogs.com/gif.latex?C_a%20%3D%20%5Csum_h%20%5Calpha_h%20%5Csum_x%20%5Cleft%5Clangle%20x%20%5Cleft%5Cvert%20D_h%20%5Cright%5Cvert%20x%20%5Cright%5Crangle%20%5Cleft%28P_%7B0x%7D%20-%20P_%7B1x%7D%20%5Cright%20%29)

where ![equation](https://latex.codecogs.com/gif.latex?D_h%20%3D%20U%5E%5Cdagger_h%20%5CSigma_h%20U_h) is the diagonal matrix and probabilities ![equation](https://latex.codecogs.com/gif.latex?P_%7Bix%7D) are obtained as a result of measurement of the following circuit

![alt text](pictures/c_vec_circ.png)

**Note:** not all ![equation](https://latex.codecogs.com/gif.latex?%5Cdpi%7B300%7D%20%5CSigma_h) have to be measured. The results for the Pauli strings differing by the replacement of Z Pauli matrices with identity matrices can be obtained from the single measurement.

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
