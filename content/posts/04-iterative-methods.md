---
title: "Iterative methods"
teaching: 10
exercises: 0
questions:
- "What are iterative methods for the solution of linear systems?"

objectives:
- ""
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
that $p_k^T r_{k-1} \ne 0$ (i.e. at each step we are guarenteed to reduce $\phi$), and 
that the search directions are linearly independent. The latter guarentees that the 
method will converge in at most $n$ steps, where $n$ is the size of the square matrix 
$A$.

It can be shown that the best set of search directions can be achieved by setting

$$
\begin{aligned}
\beta_k &= \frac{-p^T_{k-1} A r_{k-1}}{p^T_{k-1} A p_{k-1} \\
p_k &= r_{k-1} + \beta_k p_{k-1} \\
\alpha_k &= \frac{p^T_k r_{k-1}}{p^T_k A p_k}
$$

leading to the final algorithm given below (reproduced from 
[Wikipedia](https://en.wikipedia.org/wiki/Conjugate_gradient_method):
       
![Conjugate Gradient algorithm](/figs/cg_pseudocode.svg)

### Software



### Other Reading

- Golub, G. H. & Van Loan, C. F. Matrix Computations, 3rd Ed. (Johns Hopkins University 
  Press, 1996). Chapter 10 

Note: based on FD matrix in previous lession:
\vspace*{1em} 

\item {\it Solution of linear systems}

\vspace*{0.5em}
\noindent
For $N=4,8,16,32,64,128$ try the following:



\begin{enumerate}
\item Solve the linear systems using $\mathbf{U}_i=A^{-1} \mathbf{f}_i$ and
  record the time this takes on a $\log$-$\log$ graph. (Omit the case $N=128$
  and note this may take a while for $N=64$.)
\item Solve the linear systems using Gaussian elimination (corresponds to
  \textsc{Matlab}'s ``\mcode{\\}'' command). Plot the time this takes on the
  same graph.
\item ($\star$) \label {cg} Now solve the systems iteratively using \textsc{Matlab}'s
  conjugate gradients solver. How many iterations are needed for each
  problem? Explain the results for the right-hand-side $\mathbf{f}_1$. For
  the right-hand-side $\mathbf{f}_2$ what is the relationship between the
  number of iterations and $N$. How long do the computations take?
\item ($\star$) Repeat \ref{cg} using \textsc{Matlab}'s built in BICGSTAB and GMRES
  solvers.
\item ($\star$) Write a function to solve a linear system using the Jacobi method. In
  terms of $N$, how many iterations does it take to converge? (Try
  $N=4,8,16,32,64$.)
\item ($\star$) Write a function to solve a linear system using the SOR method. For
  $N=64$ and right-hand-side $\mathbf{f}_2$ determine numerically the best
  choice of the relaxation parameter to 2 decimal places and compare this
  with theory.
\end{enumerate}

\end{enumerate}
