---
title: "CIVL537 CST FEA Solver"
---

# CIVL 537 CST FEA Solver

## Overview

This is a Constant Strain Triangle (CST) Finite Element Analysis (FEA) solver
developed as part of a class project for UBC Civil 537 (Computation Mechanics).
[The live web app is on Streamlit](https://civl537-fea-solver.streamlit.app/). 

## Part 1 - Repository Set-up

+ The set-up of Docker and the Git / GitHub followed closely with the provided
  guide. 
+ Instead of Docker desktop, the docker engine and associated tools were
  installed locally to the student device running the Debian 13 (Trixie) 
  GNU/Linux OS using the corresponding [installation
  guide](https://docs.docker.com/engine/install/debian/).
+ The starter app was successfully launched locally. Deployment to Streamlit
  was also successful. 

## Part 2 - CST Element Formulation

### Area Computation

The signed area of a *parallelogram* spanned by two vectors $\langle p
\rangle=\langle p_1,p_2\rangle$ and $\langle q \rangle=\langle
q_1,q_2\rangle$ sharing a common starting point is given by: 

$$
\begin{align}
\det{\left(
\begin{bmatrix}
    p_1 & p_2 \\
    q_1 & q_2 
\end{bmatrix} \right)} 
&= \text{Area of Parallelogram Spanned}
\end{align}
$$

Where the area is positive if the shared point, the tip of $\langle p \rangle$
and the tip of $\langle q \rangle$ is ordered counterclockwise (i.e. lesser
angle from $\langle p \rangle$ to $\langle q \rangle$ is counterclockwise), and
negative if opposite [^citeFrose18]. Then, the signed area of the *triangle*
spanned by the two vectors is simply half of the above expression. in other
words:

$$
\begin{align}
\Delta &:=
\frac{1}{2}
\det{\left(
\begin{bmatrix}
    p_1 & p_2 \\
    q_1 & q_2 
\end{bmatrix} \right)} 
= \text{Area of Triangle Spanned}
\end{align}
$$

<p align="center">
<img 
    src="./readme_figs/fig_signed_area.png" 
    alt="Signed Area Illustration"
    width="500px"/>
</p>

For the CST element implementation, three points $(x_0, y_0)$, $(x_1, y_1)$ and
$(x_2, y_2)$ are given. The above formula was used on the vector from $(x_0,
y_0)$ to $(x_1, x_1)$ and the vector from $(x_0, y_0)$ to $(x_2, y_2)$, in
otherwords: 

$$
\begin{align}
\Delta &= 
\frac{1}{2}
\det\left(
\begin{bmatrix}
(x_1 - x_0) & (y_1 - y_0) \\
(x_2 - x_0) & (y_2 - y_0) \\
\end{bmatrix}
\right)
\end{align}
$$

### Shape Functions and Strain Displacement Matrix

The shape functions $[N]$ and the strain-displacement matrix
$[B]$ are derived based on the lecture notes [^citeL6Fahimi26]. For the
CST, the nodal displacement DOFs $u$ and $v$ can be written as linear
combinations of $x$ and $y$:

$$
\begin{align}
u(x, y)
    &= \alpha_1 + \alpha_2 x + \alpha_3 y \\
v(x, y)
    &= \beta_1 + \beta_2 x + \beta_3 y 
\end{align}
$$

The strains can be written as:

$$
\begin{align}
    \begin{Bmatrix}
        \epsilon_x \\
        \epsilon_y \\
        \gamma_{xy} \\
    \end{Bmatrix}
        &=
    \begin{Bmatrix}
        \frac{\partial u}{\partial x} \\
        \frac{\partial v}{\partial y} \\
        \frac{\partial u}{\partial y} + 
        \frac{\partial v}{\partial x} \\
    \end{Bmatrix}
    =
    \begin{Bmatrix}
        \alpha_2 \\
        \beta_3 \\
        \alpha_3 + \beta_2 \\
    \end{Bmatrix}
    = \text{Constant Vector}
\end{align}
$$

Then at each element's nodes $i$, $j$, $k$, we have:

$$
\begin{align}
    u_i &= \alpha_1 + \alpha_2 x_i + \alpha_3 y_i \\
    u_j &= \alpha_1 + \alpha_2 x_j + \alpha_3 y_j \\
    u_k &= \alpha_1 + \alpha_2 x_k + \alpha_3 y_k \\
    v_i &= \beta_1 + \beta_2 x_i + \beta_3 y_i \\
    v_j &= \beta_1 + \beta_2 x_j + \beta_3 y_j \\
    v_k &= \beta_1 + \beta_2 x_k + \beta_3 y_k 
\end{align}
$$

In matrix form, we can write:

$$
\begin{align}
\begin{Bmatrix}
    u_i \\ 
    u_j \\ 
    u_k \\
\end{Bmatrix}
&= 
\underbrace{
\begin{bmatrix}
1 & x_i & y_i \\
1 & x_j & y_j \\
1 & x_k & y_k \\
\end{bmatrix}}_{[A]}
\underbrace{
\begin{Bmatrix}
    \alpha_1 \\ 
    \alpha_2 \\ 
    \alpha_3 \\
\end{Bmatrix}}_{\begin{Bmatrix}\alpha\end{Bmatrix}} \\
\begin{Bmatrix}
    v_i \\ 
    v_j \\ 
    v_k \\
\end{Bmatrix}
&= 
\underbrace{
\begin{bmatrix}
1 & x_i & y_i \\
1 & x_j & y_j \\
1 & x_k & y_k \\
\end{bmatrix}}_{[A]}
\underbrace{
\begin{Bmatrix}
    \beta_1 \\ 
    \beta_2 \\ 
    \beta_3 \\
\end{Bmatrix}}_{\begin{Bmatrix}\beta\end{Bmatrix}}
\end{align}
$$

Solving for the coefficients $\alpha_n$ and $\beta_n$, we first invert the
matrix $[A]$:

$$
\begin{align}
    [A]^{-1}
    &= \frac{1}{\det[A]}
    \begin{bmatrix}
        x_j y_k - x_k y_j & x_k y_i - x_i y_k & x_i y_j - x_j y_i \\
        y_j - y_k & y_k - y_i & y_i - y_j \\
        x_k - x_j & x_i - x_k & x_j - x_i \\
    \end{bmatrix} 
    = \frac{1}{2\Delta}
    \begin{bmatrix}
        a_i & a_j & a_k \\
        b_i & b_j & b_k \\
        c_i & c_j & c_k \\
    \end{bmatrix}
\end{align}
$$

Where the signed area of the triangle $\Delta$ was previously defined. Then
solving the coefficients, we have: 

$$
\begin{align}
\begin{Bmatrix}
    \alpha_1 \\ 
    \alpha_2 \\ 
    \alpha_3 \\
\end{Bmatrix}
&= 
\frac{1}{2\Delta}
\begin{bmatrix}
    a_i & a_j & a_k \\
    b_i & b_j & b_k \\
    c_i & c_j & c_k \\
\end{bmatrix}
\begin{Bmatrix}
    u_i \\ 
    u_j \\ 
    u_k \\
\end{Bmatrix} \\
\begin{Bmatrix}
    \beta_1 \\ 
    \beta_2 \\ 
    \beta_3 \\
\end{Bmatrix}
&=
\frac{1}{2\Delta}
\begin{bmatrix}
    a_i & a_j & a_k \\
    b_i & b_j & b_k \\
    c_i & c_j & c_k \\
\end{bmatrix}
\begin{Bmatrix}
    v_i \\ 
    v_j \\ 
    v_k \\
\end{Bmatrix} \\
\end{align}
$$

Then substituting these coefficients back into 
$u(x, y) = \alpha_1 + \alpha_2 x + \alpha_3 y$ 
and $v(x, y) = \beta_1 + \beta_2 x + \beta_3 y$, then gathering terms by $u_n$
and $v_n$, we have:

$$
\begin{align}
u(x, y) 
    &=
    \sum_{n=i,j,k}
    \underbrace{
    \frac{1}{2\Delta}
    (a_n + b_n x + c_n)}_{N_n} u_n 
    = \sum_{n=i,j,k} N_n u_n \\
v(x, y) 
    &=
    \sum_{n=i,j,k}
    \underbrace{
    \frac{1}{2\Delta}
    (a_n + b_n x + c_n)}_{N_n} v_n
    = \sum_{n=i,j,k} N_n v_n \\
N_n
    &= \frac{1}{2\Delta} (a_n + b_n + c_n)
\end{align}
$$

In matrix form, the ***shape functions*** $[N]$ can be written: 

$$
\begin{align}
\begin{Bmatrix}
    u(x, y) \\
    v(x, y) \\
\end{Bmatrix}
&=
\underbrace{
\begin{bmatrix}
    N_i & 0 & N_j & 0 & N_k & 0 \\
    0 & N_i & 0 & N_j & 0 & N_k \\
\end{bmatrix}
}_{[N]}
\underbrace{
\begin{Bmatrix}
    u_i \\
    v_i \\
    u_j \\
    v_j \\
    u_k \\
    v_k \\
\end{Bmatrix}
}_{\begin{Bmatrix}s\end{Bmatrix}}
= [N]\begin{Bmatrix}s\end{Bmatrix}
\end{align}
$$

Then, we can compute the ***strain-displacement matrix*** $[B]$ using
the definition of the strains: 

$$
\begin{align}
    \begin{Bmatrix}
        \epsilon_x \\
        \epsilon_y \\
        \gamma_{xy} \\
    \end{Bmatrix}
    &=
    \begin{Bmatrix}
        \frac{\partial u}{\partial x} \\
        \frac{\partial v}{\partial y} \\
        \frac{\partial u}{\partial y} + 
        \frac{\partial v}{\partial x} \\
    \end{Bmatrix}
    =
    \begin{Bmatrix}
        \alpha_2 \\
        \beta_3 \\
        \alpha_3 + \beta_2 \\
    \end{Bmatrix} \\
    \begin{Bmatrix}
        \epsilon_x \\
        \epsilon_y \\
        \gamma_{xy} \\
    \end{Bmatrix}
    &= \frac{1}{2\Delta}
    \begin{Bmatrix}
    b_i u_i + b_j u_j + b_k u_k \\
    c_i v_i + c_j v_j + c_k v_j \\
    c_i u_i + c_j u_j + c_k u_k + b_i v_i + b_j v_j + b_k v_k \\
    \end{Bmatrix} \\
    \begin{Bmatrix}
        \epsilon_x \\
        \epsilon_y \\
        \gamma_{xy} \\
    \end{Bmatrix}
    &= 
    \underbrace{
    \frac{1}{2\Delta}
    \begin{bmatrix}
    b_i & 0 & b_j & 0 & b_k & 0 \\
    0 & c_i & 0 & c_j & 0 & c_k \\
    c_i & b_i & c_j & b_j & c_k & b_k \\
    \end{bmatrix}}_{[B]}
    \begin{Bmatrix}
    u_i \\
    v_i \\
    u_j \\
    v_j \\
    u_k \\
    v_k \\
    \end{Bmatrix} = [B]\begin{Bmatrix}s\end{Bmatrix} 
\end{align}
$$

The above results for $[B]$ is used for the element implementation of
the CSTs.

### Constitutive Matrix

The strains and stresses involved in the strain and stress tensors of 3D solids
are illustrated as follows:

<p align="center">
<img 
    src="./readme_figs/fig_tensor_dirs.png" 
    alt="Strains and Stresses of a 3D Solid"
    width=1000px>
</p>

For an isotropic material in equilibrium, it requires $\tau_{ij}=\tau_{ji}$ for
the shear stresses and $\epsilon_{ij} = \epsilon_{ji}$ for the shear strains.
Further, the shear strains are often re-written in terms of ***engineering
shear strain*** as follows [^citeC2_2Bowers25]:

$$
\gamma_{ij}=\gamma_{ji} = 2\epsilon_{ij} = 2\epsilon_{ji}
$$

<p align="center">
<img 
    src="./readme_figs/fig_shear_strain.png" 
    alt="Shear Strain and Engineering Shear Strain"
    width=500px/>
</p>

The constitutive relations for an isotropic material are[^citeC3_2Bowers25]:

$$
\begin{align}
\epsilon_{ii} 
&= \frac{\sigma_{ii}}{E} 
    - \nu \frac{\sigma_{jj}}{E} 
    - \nu \frac{\sigma_{kk}}{E} \\
\gamma_{ij}
&= \frac{\tau_{ij}}{G}
\end{align}
$$

Where $E$ is the ***Elastic modulus***, $\nu$ is the ***Poisson's ratio***, and
$G$ is the ***shear modulus***, which for an isotropic material is:

$$
\begin{align}
G &= \frac{E}{2(1+\nu)}
\end{align}
$$

In matrix form, we have: 

$$
\begin{align}
\begin{Bmatrix}
    \epsilon_{xx} \\
    \epsilon_{yy} \\
    \epsilon_{zz} \\
    \gamma_{yz} \\
    \gamma_{xz} \\
    \gamma_{xy} \\
\end{Bmatrix}
&=
\frac{1}{E}
\begin{bmatrix}
    1 & -\nu & -\nu & 0 & 0 & 0 \\
    -\nu & 1 & -\nu & 0 & 0 & 0 \\
    -\nu & -\nu & 1 & 0 & 0 & 0 \\
    0 & 0 & 0 & 2(1+\nu) & 0 & 0 \\
    0 & 0 & 0 & 0 & 2(1+\nu) & 0 \\
    0 & 0 & 0 & 0 & 0 & 2(1+\nu) \\
\end{bmatrix}
\begin{Bmatrix}
    \sigma_{xx} \\
    \sigma_{yy} \\
    \sigma_{zz} \\
    \tau_{yz} \\
    \tau_{xz} \\
    \tau_{xy} \\
\end{Bmatrix}
\end{align}
$$

Taking the inverse, we will obtain the general constitutive matrix for an
isotropic 3D material: 

$$
\begin{align}
\begin{Bmatrix}
    \sigma_{xx} \\
    \sigma_{yy} \\
    \sigma_{zz} \\
    \tau_{yz} \\
    \tau_{xz} \\
    \tau_{xy} \\
\end{Bmatrix}
&=
\underbrace{
\frac{E}{(1+\nu)(1-2\nu)}
\begin{bmatrix}
    1-\nu & \nu & \nu & 0 & 0 & 0 \\
    \nu & 1-\nu & \nu & 0 & 0 & 0 \\
    \nu & \nu & 1-\nu & 0 & 0 & 0 \\
    0 & 0 & 0 & \frac{1-2\nu}{2} & 0 & 0 \\
    0 & 0 & 0 & 0 & \frac{1-2\nu}{2} & 0 \\
    0 & 0 & 0 & 0 & 0 & \frac{1-2\nu}{2} \\
\end{bmatrix}}_{[D]}
\underbrace{
\begin{Bmatrix}
    \epsilon_{xx} \\
    \epsilon_{yy} \\
    \epsilon_{zz} \\
    \gamma_{yz} \\
    \gamma_{xz} \\
    \gamma_{xy} \\
\end{Bmatrix}}_{\begin{Bmatrix}\epsilon\end{Bmatrix}}
= [D]\begin{Bmatrix}\epsilon\end{Bmatrix}
\end{align}
$$

By setting $\sigma_{zz}=\tau_{yz}=\tau_{xz}=0$, we obtain the ***Plane Stress
Constitutive Matrix*** as provided in the lecture
[^citeL6Fahimi26] [^citeC3_2Bowers25]:

$$
\begin{align}
[D]
&=
\frac{E}{1-\nu^2}
\begin{bmatrix}
1 & \nu & 0 \\
\nu & 1 & 0 \\
0 & 0 & \frac{1-\nu}{2} \\
\end{bmatrix}
\end{align}
$$

On the other hand, by setting $\epsilon_{zz}=\epsilon_{yz}=\epsilon_{xz}=0$,
the ***Plane Strain Constitutive Matrix*** is obtained:

$$
\begin{align}
[D]
&=
\frac{E(1-\nu)}{(1+\nu)(1-2\nu)}
\begin{bmatrix}
1 & \frac{\nu}{1-\nu} & 0 \\
\frac{\nu}{1-\nu} & 1 & 0 \\
0 & 0 & \frac{1-2\nu}{2(1-\nu)} \\
\end{bmatrix}
\end{align}
$$

The above matrices were implemented accordingly. 

### Stiffness Matrix

Based on the principle of virtual work formulation, the element stiffness
matrix $[k_e]$ is:

$$
\begin{align}
[k_e] 
&= \int_{V_e} [B]^T [D] [B] dV \\
[k_e] 
&= \iiint [B]^T [D] [B] dx dy dz
\end{align}
$$

For the CST element, based on the earlier derived $[B]$ and $[D]$, we see that
they do not depend on the variables of integration (only pre-defined nodal
coordinates and material properties). If the thickness of the elemnt is $h$,
then the expression simplifies to:

$$
\begin{align}
[k_e] 
&= [B]^T [D] [B] h \iint dx dy \\
[k_e] 
&= [B]^T [D] [B] h \Delta
\end{align}
$$

Which provides the element stiffness matrix for implementation. 

### Hand Verification Example

The strain-displacement matrix $[B]$ is verified using hand calculation with
the following points: 

$$
\begin{align}
(x_i, y_i) &= (0, 0) \\
(x_j, y_j) &= (1, 0) \\
(x_k, y_k) &= (0, 1) 
\end{align}
$$

Computing coefficients: 

$$
\begin{align}
a_i 
    &= x_j y_k - x_k y_j = 1 \times 1 - 0 \times 0 = 1 \\
a_j 
    &= x_k y_i - x_i y_k = 0 \times 0 - 0 \times 1 = 0 \\
a_k 
    &= x_i y_j - x_j y_i = 0 \times 0 - 1 \times 0 = 0 \\
b_i 
    &= y_j - y_k = 0 - 1 = -1 \\
b_j 
    &= y_k - y_i = 1 - 0 = 1 \\
b_k 
    &= y_i - y_j = 0 - 0 \\
c_i 
    &= x_k - x_j = 0 - 1 = -1 \\
c_j 
    &= x_i - x_k = 0 - 0 = 0 \\
c_k 
    &= x_j - x_i = 1 - 0 = 1
\end{align}
$$

Compute area:

$$
\begin{align}
\Delta
    &= \frac{1}{2}\det\left(
    \begin{bmatrix}
    x_j - x_i & y_j - y_i \\
    x_k - x_i & y_k - y_i \\
    \end{bmatrix}
    \right) \\
\Delta
    &= \frac{1}{2}\det\left(
    \begin{bmatrix}
    1 & 0 \\
    0 & 1 \\
    \end{bmatrix}
    \right) = \frac{1}{2} \\
2\Delta
    &= 1
\end{align}
$$

Computing shape functions: 

$$
\begin{align}
N_n 
    &= \frac{1}{2\Delta}(a_n + b_n x + c_n y) \\
\implies
N_i 
    &= 1 - x - y \\
N_j 
    &= x \\
N_k 
    &= y 
\end{align}
$$

Check shape functions:

$$
\begin{align}
\sum_{n=i,j,k} N_n 
    &= 1 - x - y + x + y = 1 \implies \text{ok} \\
\end{align}
$$

$$
\begin{align}
    N_i(0,0) &= 1 + 0 + 0 = 1 \implies \text{ok} \\
    N_j(0,0) &= 0 \implies \text{ok} \\
    N_k(0,0) &= 0 \implies \text{ok} 
\end{align}
$$

$$
\begin{align}
    N_i(1,0) &= 1 - 1 + 0 = 0 \implies \text{ok} \\
    N_j(1,0) &= 1 \implies \text{ok} \\
    N_k(1,0) &= 0 \implies \text{ok} 
\end{align}
$$

$$
\begin{align}
    N_i(0,1) &= 1 - 0 + 1 = 0 \implies \text{ok} \\
    N_j(0,1) &= 0 \implies \text{ok} \\
    N_k(0,1) &= 1 \implies \text{ok} 
\end{align}
$$

Differentiate shape functions:

$$
\begin{align}
\begin{matrix}
\frac{\partial N_i}{\partial x}
    = -1 
& \frac{\partial N_j}{\partial x}
    = 1 
& \frac{\partial N_k}{\partial x}
    = 0 \\
\frac{\partial N_i}{\partial y}
    = -1 
& \frac{\partial N_j}{\partial y}
    = 0 
& \frac{\partial N_k}{\partial y}
    = 1 \\
\end{matrix}
\end{align}
$$

Compute strain vector and strain-displacement matrix:

$$
\begin{align}
\begin{Bmatrix}
\epsilon_x \\
\epsilon_y \\
\gamma_{xy}
\end{Bmatrix}
&=
\begin{Bmatrix}
\frac{\partial u}{\partial x} \\
\frac{\partial v}{\partial y}  \\
\frac{\partial u}{\partial y} 
    +\frac{\partial v}{\partial x} \\
\end{Bmatrix} =
\begin{Bmatrix} 
\frac{\partial N_n}{\partial x} u_n \\
\frac{\partial N_n}{\partial y} v_n \\
\frac{\partial N_n}{\partial y} u_n 
    +\frac{\partial N_n}{\partial x} v_n \\
\end{Bmatrix} \\
\begin{Bmatrix}
\epsilon_x \\
\epsilon_y \\
\gamma_{xy}
\end{Bmatrix}
&=
\begin{Bmatrix}
-u_i + u_j \\
-v_i + v_k \\
-u_i - v_i + v_j + u_k \\
\end{Bmatrix} \\
\begin{Bmatrix}
\epsilon_x \\
\epsilon_y \\
\gamma_{xy}
\end{Bmatrix}
&=
\underbrace{
\begin{bmatrix}
-1 & 0 & 1 & 0 & 0 & 0 \\
0 & -1 & 0 & 0 & 0 & 1 \\
-1 & -1 & 0 & 1 & 1 & 0 \\
\end{bmatrix}}_{[B]}
\underbrace{
\begin{Bmatrix}
u_i \\
v_i \\
u_j \\
v_j \\
u_k \\
v_k \\
\end{Bmatrix}}_{\begin{Bmatrix} s \end{Bmatrix}}
\end{align}
$$

Compare with computation results:

```python
>>> from src.elements import compute_D, compute_B, compute_area, compute_k
>>> import numpy as np
>>> coords = np.array([[0, 0], [1, 0], [0, 1]])
>>> compute_B(coords)
array([[-1.,  0.,  1.,  0.,  0.,  0.],
       [ 0., -1.,  0.,  0.,  0.,  1.],
       [-1., -1.,  0.,  1.,  1.,  0.]])
```

We see that the hand calculation and the implemented function aligns. 

## Part 3 - Meshing

## Part 4 - Global Assembly

### Global Stiffness Matrix

There are two DOFs per global node (before applying B.C.s). A simple global DOF
indexing scheme is adopted, where for each global node number $n$, the
corresponding DOF indices are:

+ For horizontal displacement: $2n$
+ For the vertical displacement: $2n+1$

The global node indices of each element is earlier stored in `elements`. If an
element has vertices at nodes $n_i$, $n_j$ and $n_k$ in sequence, the
element-to-global DOF indices will be mapped as:

| Element DOF Index $\mathrm{el}$ | Global DOF Index $\mathrm{gl}$ | DOF Description             |
| -------                         | --------                       | ------                      |
| $0$                             | $2n_i$                         | Horizontal DOF at i-th Node |
| $1$                             | $2n_i+1$                       | Vertical DOF at i-th Node   |
| $2$                             | $2n_j$                         | Horizontal DOF at j-th Node |
| $3$                             | $2n_j+1$                       | Vertical DOF at j-th Node   |
| $4$                             | $2n_k$                         | Horizontal DOF at k-th Node |
| $5$                             | $2n_k+1$                       | Vertical DOF at k-th Node   |

So to assemble the global stiffness matrix $[K_g]$, the following procedure is
applied *for each element*:

1. Get the coordinates for the element's nodes $n_i$, $n_j$, $n_k$.
2. Construct the element stiffness matrix $[k_e]$ using the coordinates at
   nodes $n_i$, $n_j$, $n_k$. 
3. Iterate over the element stiffness matrix rows and columns. Using the above
   index mapping of $\mathrm{gl}(\mathrm{el})$, add the entries in the element
   stiffness matrix to the global stiffness matrix (where subscript $r$, $c$
   means DOF indices for the row and column respectively).
   $$ 
   [K_g](\mathrm{gl}(\mathrm{el}_r),\mathrm{gl}(\mathrm{el}_c)) {+=} [k_e](\mathrm{el}_r, \mathrm{el}_c)
   $$

# References

[^citeFrose18]: R. Froese and B. Wetton, "Notes for Math 152: Linear Systems",
    2018.  [Online]. Available
    [https://personal.math.ubc.ca/~karu/m152/notes.pdf](https://personal.math.ubc.ca/~karu/m152/notes.pdf) 
[^citeL6Fahimi26]: S. Fahimi. (2026). UBC CIVL 537 Computation Mechanics: "06 -
    Continuum Finite Elements"
[^citeC2_2Bowers25]: A. F. Bowers, "Mathematical description of shape changes
    in solids" in *Applied Mechanics of Solids*, 2025, ch. 2.2. [Online].
    Available:
    [https://solidmechanics.org/Text/Chapter2_2/Chapter2_2.php](https://solidmechanics.org/Text/Chapter2_2/Chapter2_2.php)
[^citeC3_2Bowers25]: A. F. Bowers, "Linear elastic material behavior" in *Applied Mechanics of Solids*, 2025, ch. 3.2. [Online].
    Available:
    [https://solidmechanics.org/Text/Chapter3_2/Chapter3_2.php](https://solidmechanics.org/Text/Chapter3_2/Chapter3_2.php)
