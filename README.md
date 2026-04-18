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
  installed locally to the student device running the Debian 13 (Bookworm) 
  GNU/Linux OS using the corresponding [installation
  guide](https://docs.docker.com/engine/install/debian/).
+ The starter app was successfully launched locally. Deployment to Streamlit
  was also successful. 

## Part 2 - CST Element Formulation

### Area Computation

The signed area of a *parallelogram* spanned by two vectors
$\mathbf{p}=<p_1,\,p_2>$ and $\mathbf{q}=<q_1,\,q_2>$ sharing a common starting
point is given by: 

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

Where the area is positive if the shared point, the tip of $\mathbf{p}$ and the
tip of $\mathbf{q}$ is ordered counterclockwise (i.e. lesser angle from
$\mathbf{p}$ to $\mathbf{q}$ is counterclockwise), and negative if opposite
[^citeFrose18]. Then, the signed area of the *triangle* spanned by the two vectors
is simply half of the above expression. in other words:

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

## Shape Functions and Strain Displacement Matrix

The shape functions $\mathbf{N}$ and the strain-displacement matrix
$\mathbf{B}$ are derived based on the lecture notes [^citeL6Fahimi26]. For the
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
\end{bmatrix}}_{\mathbf{A}}
\underbrace{
\begin{Bmatrix}
    \alpha_1 \\ 
    \alpha_2 \\ 
    \alpha_3 \\
\end{Bmatrix}}_{\mathbf{\alpha}} \\
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
\end{bmatrix}}_{\mathbf{A}}
\underbrace{
\begin{Bmatrix}
    \beta_1 \\ 
    \beta_2 \\ 
    \beta_3 \\
\end{Bmatrix}}_{\mathbf{\beta}}
\end{align}
$$

Solving for the coefficients $\alpha_n$ and $\beta_n$, we first invert the
matrix $\mathbf{A}$:

$$
\begin{align}
    \mathbf{A}^{-1}
    &= \frac{1}{\det\mathbf{A}}
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
    (a_n + b_n x + c_n)}_{N_n} u_n \\
v(x, y) 
    &=
    \sum_{n=i,j,k}
    \underbrace{
    \frac{1}{2\Delta}
    (a_n + b_n x + c_n)}_{N_n} v_n
\end{align}
$$


# References

[^citeFrose18]: R. Froese and B. Wetton, "Notes for Math 152: Linear Systems",
    2018.  [Online]. Available
    [https://personal.math.ubc.ca/~karu/m152/notes.pdf](https://personal.math.ubc.ca/~karu/m152/notes.pdf) 
[^citeL6Fahimi26]: S. Fahimi, "06 - Continuum Finite Elements", 2026.
