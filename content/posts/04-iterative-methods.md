---
title: "Iterative methods"
teaching: 10
exercises: 0
questions:
- "What are iterative methods for the solution of linear systems?"

objectives:
- ""
katex: true
markup: "mmark"
---

Previously we have discussed *direct* linear algebra solvers based on decompositions of 
the original matrix $A$. The amount of computational effort required to achieve these 
decomposisions is $\mathcal{O}(n^3)$, where $n$ is the number of rows of a square 
matrix. They are therefore unsuitable for the large, sparse systems of equations that 
are typically encountered in scientific applications. An alternate class of linear 
algebra solvers are the *iterative* methods, which produce a series of *approximate* 
solutions $x_k$ to the $A x = b$ problem. The performance of each algorithm is then 
based on how quickly, or how many iterations $k$ are required, for the solution $x_k$ to 
converge to within a set tolerance of the true solution $x$.

## Jacobi Method

The Jacobi method is the simplest of the iterative methods, and relies on the fact that 
the matrix is *diagonally dominant*. Starting from the problem definition:

$$
A\mathbf{x} = \mathbf{b}
$$

we decompose $A$ in to $A = L + D + U$, where $L$ is lower triagular, $D$ is diagonal, 
$U$ is upper triangular. 

$$
A\mathbf{x} = L\mathbf{x} + D\mathbf{x} + U\mathbf{x} =  \mathbf{b}
$$

We then assume that we have an initial guess at the solution $\mathbf{x}^0$, and try to 
find a new estimate $\mathbf{x}^1$. Assuming that the diagonal $D$ dominates over $L$ 
and $U$, a sensible choice would be to insert $x^0$ and the unknown $x^1$ into the 
equation like so:

$$
L\mathbf{x}^0 + D\mathbf{x}^1 + U\mathbf{x}^0 =  \mathbf{b}
$$

we can rearrange to get an equation for $x^1$. This is easily solved as we can take the 
inverse of the diagonal matrix by simply inverting each diagonal element individually:

$$
D\mathbf{x}_1 =  \mathbf{b} - (L+U)\mathbf{x}_0
$$

Thus we end up with the general Jacobi iteration:

$$
\mathbf{x}_{k+1} =  D^{-1}(\mathbf{b} - (L+U)\mathbf{x}_k)
$$

## Relaxation methods

The Jacobi method is an example of a relaxation method, where the matrix $A$ is split 
into a dominant part $M$ (which is easy to solve), and the remainder $N$. That is, $A = 
M - N$

$$M\mathbf{x}_{k+1} = N\mathbf{x}_k + \mathbf{b}$$
$$\mathbf{x}_{k+1} = M^{-1}N\mathbf{x}_k + M^{-1}\mathbf{b}$$

For the Jacobi method $M = D$ and $N = -(L + U)$. Other relaxation methods include 
Gauss-Seidel, where $M = (D + L)$ and $N = -U$, and successive over-relaxation (SOR), 
where $M = D + \omega L$ and $N = (1 - \omega) D - \omega U$, where $\omega$ is the 
*relaxation* parameter.

For any relaxation method to converge we need $\rho(M^{-1}N) < 1$, where $\rho()$ is the 
*spectral radius* of $M^{-1} N$, which is defined as the largest eigenvalue $\lambda$ of 
a a given matrix $G$:

$$
\rho(G) = \max{|\lambda|: \lambda \in \lambda(G)}
$$

For the SOR method, the relaxation parameter $\omega$ is generally chosen to minimise 
$\rho(M^{-1}N)$, so that the speed of convergence is maximised. In some cases this 
optimal $\omega$ is known, for example for finite difference discretisation of the 
[Poisson equation](https://www.sciencedirect.com/science/article/pii/S0893965908001523).
However, in many cases sophisticated eigenvalue analysis is required to determine the 
optimal $\omega$. 

### Problems

This exercise involves the manipulation and solution of the linear system resulting from 
the finite difference solution to Poisson's equation in two dimensions. Let $A$ be a 
sparse symmetric positive definite matrix of dimension $(N-1)^2 \times (N-1)^2$ created 
using `scipy.sparse` (for a given $N$) by the function
`buildA` as follows:
```python
import numpy as np
import scipy.sparse as sp

def buildA(N):
  dx = 1 / N
  nvar = (N - 1)^2;
  e1 = np.ones((nvar, 1));
  e2 = e1
  e2[1:N-1:nvar] = 0
  e3 = e1
  e3[N-1:N-1:nvar] = 0
  A = sp.spdiags(
        np.vstack((-e1, 4*e1, -e1)),
        -(N-1):N-1:N-1, nvar, nvar
      ) +
      sp.spdiags(
        np.vstack((-e3, -e2),
        -1:2:1 , nvar, nvar
      )
  A= A / dx^2;
```

and let $\mathbf{f}_1$ and $\mathbf{f}_2$ be the vectors defined in
`buildf1` and `buildf2`

```python
def buildf1(N):
  x = 0:1/N:1
  y = x
  f = np.dot(np.sin(pi*x), np.sin(pi*y))
  return f[2:N,2:N].reshape(-1,1)
```

```python
def buildf2(N):
  x = 0:1/N:1
  y = x
  f = np.dot(np.max(x,1-x), np.max(y,1-y))
  return f[2:N,2:N].reshape(-1, 1)
```

We will consider manipulation of the matrix $A$ and solution of the linear
systems $A\mathbf{U}_i=\mathbf{f}_i$. The solution to this linear system
corresponds to a finite difference solution to Poisson's equation $-\nabla^2 u
= f$ on the unit square with zero Dirichlet boundary conditions where $f$ is
either $\sin(\pi x) \sin (\pi y)$ or $\max(x,1-x) \max(y,1-y)$. PDEs of this type occur 
(usually with some additional reaction and or convection terms) very frequently
in mathematical modelling of physiological processes, and even in image
analysis. 

1. Write a function to solve a linear system using the Jacobi method. In
  terms of $N$, how many iterations does it take to converge? (Try
  $N=4,8,16,32,64$.)
2. Write a function to solve a linear system using the SOR method. For
  $N=64$ and right-hand-side $\mathbf{f}_2$ determine numerically the best
  choice of the relaxation parameter to 2 decimal places and compare this
  with theory.

## Conjugate Gradient Method

One of the most important classes of iterative methods are the *Krylov subspace 
methods*, which include:
- *Congugate Gradient (CG)*: for symmetrix positive definite matrices
- *Biconjugate Gradient Stabilized (BiCGSTAB)*: for general square matrices
- *Generalized Minimal Residual (GMRES)*: for general square matrices

Below we will give a brief summary of the CG method, for more details you can consult 
the text by Golub and Van Loan (Chapter 10).

The CG method is based on minimising the function

$$
\phi(x) = \frac{1}{2}x^T A x - x^T b
$$

If we set $x$ to the solution of $Ax =b$, that is $x = A^{-1} b$, then the value of 
$\phi(x)$ is at its minimum $\phi(A^{-1} b) = -b^T A^{-1} b / 2$, showing that solving 
$Ax = b$ and minimising $\phi$ are equivalent.

At each iteration $k$ of CG we are concerned with the *residual*, defined as $r_k = b - 
A x_k$. If the residual is nonzero, then at each step we wish to find a positive 
$\alpha$ such that $\phi(x_k + \alpha p_k) < \phi(x_k)$, where $p_k$ is the *search 
direction* at each $k$. For the classical stepest descent optimisation algorithm the 
search direction would be the residual $p_k = r_k$, however, steapest descent can suffer 
from convergence problems, so instead we aim to find a set of search directions $p_k$ so 
that $p_k^T r\_{k-1} \ne 0$ (i.e. at each step we are guarenteed to reduce $\phi$), and 
that the search directions are linearly independent. The latter guarentees that the 
method will converge in at most $n$ steps, where $n$ is the size of the square matrix 
$A$.

It can be shown that the best set of search directions can be achieved by setting

$$
\begin{aligned}
\beta_k &= \frac{-p^T_{k-1} A r_{k-1}}{p^T_{k-1} A p_{k-1}} \\
p_k &= r_{k-1} + \beta_k p_{k-1} \\
\alpha_k &= \frac{p^T_k r_{k-1}}{p^T_k A p_k}
\end{aligned}
$$

Directly using the above equations in an iterative algorithm results in the standard CG 
algorithm. A more efficient algorithm can be derived from this by computing the 
residuals recursivly via $r_k = r\_{k-1} - \alpha_k A p_k$, leading to the final 
algorithm given below (reproduced from 
[Wikipedia](https://en.wikipedia.org/wiki/Conjugate_gradient_method)):

![Conjugate Gradient algorithm](/figs/cg_pseudocode.svg)

### Preconditioning

The CG method works well (i.e. converges quickly) if the *condition number* of the 
matrix $A$ is low. The condition number of a matrix gives a measure of how much the 
solution $x$ changes in response to a small change in the input $b$, and is a property 
of the matrix $A$ itself, so can vary from problem to problem. In order to keep the 
number of iterations small for iterative solvers, it is therefore often neccessary to 
use a *preconditioner*, which is a method of transforming what might be a difficult 
problem with a poorly conditioned $A$, into a well conditioned problem that is easy to 
solve.

Consider the case of precoditioning for the CG methods, we start from the standard 
problem $A x = b$, and we wish to solve an *equivilent* transformed problem given by

$$
\tilde{A} \tilde{x} = \tilde{b}
$$

where $\tilde{A} = C^{-1} A C^{-1}$, $\tilde{x} = Cx$, $\tilde{b} = C^{-1}$, and $C$ is 
a symmetric positive matrix.

We then simply apply the standard CG method as given above to this transformed problem. 
This leads to an algorithm which is then simplified by instead computing the transformed 
quantities $\tilde{p}_k = C p_k$, $\tilde{x}_k = C x_k$, and $\tilde{r}_k = C^{-1} r_k$. 
Finally we define a matrix $M = C^2$, which is known as the *preconditioner*, leading to 
the final precoditioned CG algorithm given below (reproduced and edited from 
[Wikipedia](https://en.wikipedia.org/wiki/Conjugate_gradient_method)):

$\mathbf{r}\_0 := \mathbf{b} - \mathbf{A x}\_0$\\
$\mathbf{z}\_0 := \mathbf{M}^{-1} \mathbf{r}\_0$\\
$\mathbf{p}\_0 := \mathbf{z}\_0$\\
$k := 0 \, $\\
**repeat until $|| \mathbf{r}_k ||_2 < \epsilon ||\mathbf{b}||_2$**\\
&nbsp;&nbsp; $\alpha\_k := \frac{\mathbf{r}\_k^T \mathbf{z}\_k}{ \mathbf{p}\_k^T 
\mathbf{A p}\_k }$\\
&nbsp;&nbsp; $\mathbf{x}\_{k+1} := \mathbf{x}\_k + \alpha\_k \mathbf{p}\_k$ \\
&nbsp;&nbsp; $\mathbf{r}\_{k+1} := \mathbf{r}\_k - \alpha_k \mathbf{A p}\_k$ \\
&nbsp;&nbsp; **if** $r\_{k+1}$ is sufficiently small then exit loop **end if** \\
&nbsp;&nbsp; $\mathbf{z}\_{k+1} := \mathbf{M}^{-1} \mathbf{r}\_{k+1}$\\
&nbsp;&nbsp; $\beta\_k := \frac{\mathbf{r}\_{k+1}^T \mathbf{z}\_{k+1}}{\mathbf{r}\_k^T 
\mathbf{z}\_k}$\\
&nbsp;&nbsp; $\mathbf{p}\_{k+1} := \mathbf{z}\_{k+1} + \beta_k \mathbf{p}\_k$\\
&nbsp;&nbsp; $k := k + 1 \, $\\
**end repeat**

The key point to note here is that the preconditioner is used by inverting $M$, so this 
matrix must be "easy" to solve in some fashion, and also result in a transformed problem 
with better conditioning.

**Termination**: The CG algorithm is normally run untill convergence to a given 
tolerance which is based on the norm of the input vector $b$. In the algorithm above we 
iterate until the residual norm is less than some fraction (set by the user) of the norm 
of $b$.

What preconditioner to choose for a given problem is often highly problem-specific, but 
some useful general purpose preconditioners exist, such as the *incomplete Cholesky 
preconditioner* for preconditioned CG (see Chapter 10.3.2 of the Golub & Van Loan text 
given below). Chapter 3 of the [Barrett et al. 
text](https://www.netlib.org/templates/templates.pdf), also cited below, contains 
descriptions of a few more commonly used preconditioners.

### Software

Once again the best resource for Python is the [`scipi.sparse.linalg` 
documentation](https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html).

Iterative methods for linear equation systems in `scipy.sparse.linalg`:

- [BIConjugate Gradient iteration 
  (BiCG)](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.bicg.html#scipy.sparse.linalg.bicg)
- [BIConjugate Gradient STABilized iteration 
  (BiCGSTAB)](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.bicgstab.html#scipy.sparse.linalg.bicgstab)
- [Conjugate Gradient iteration 
  (CG)](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.cg.html#scipy.sparse.linalg.cg)
- [Conjugate Gradient Squared iteration 
  (CGS)](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.cgs.html#scipy.sparse.linalg.cgs)
- [Generalized Minimal RESidual iteration 
  (GMRES)](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.gmres.html#scipy.sparse.linalg.gmres)
- 
  [LGMRES](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.lgmres.html#scipy.sparse.linalg.lgmres)
- [MINimum RESidual iteration 
  (MINRES)](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.minres.html#scipy.sparse.linalg.minres)
- [Quasi-Minimal Residual iteration 
  (QMR)](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.qmr.html#scipy.sparse.linalg.qmr)
- 
  [GCROT(m,k)](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.gcrotmk.html#scipy.sparse.linalg.gcrotmk)

`scipy.sparse.linalg` also contains two iterative solvers for least-squares problems, 
[`lsqr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.lsqr.html#scipy.sparse.linalg.lsqr) 
and 
[`lsmr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.lsmr.html#scipy.sparse.linalg.lsmr)

### Other Reading

- Golub, G. H. & Van Loan, C. F. Matrix Computations, 3rd Ed. (Johns Hopkins University 
  Press, 1996). Chapter 10 
- Barrett, R., Berry, M., Chan, T. F., Demmel, J., Donato, J., Dongarra, J., ... & Van 
  der Vorst, H. (1994). Templates for the solution of linear systems: building blocks 
  for iterative methods. Society for Industrial and Applied Mathematics.

### Problems

Note: based on the FD matrix in previous exercise:

For $N=4,8,16,32,64,128$ try the following:
1. Solve the linear systems using $\mathbf{U}_i=A^{-1} \mathbf{f}_i$ (see 
   [`scipy.linalg.inv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.inv.html) 
   and record the time this takes on a $\log$-$\log$ graph. (Omit the case $N=128$
  and note this may take a while for $N=64$.)
2. Solve the linear systems using $LU$ and $Cholesky$ decomposition. Plot the time this 
   takes on the same graph.
3. Now solve the systems iteratively using a conjugate gradients solver (you can use the 
   one in `scipy.linalg.sparse`, or you can code up your own). How many iterations are 
   needed for each problem? Explain the results for the right-hand-side $\mathbf{f}_1$. 
   For the right-hand-side $\mathbf{f}_2$ what is the relationship between the number of 
   iterations and $N$. How long do the computations take?
4. Repeat using the `scipy.sparse.linalg` BICGSTAB and GMRES solvers.
