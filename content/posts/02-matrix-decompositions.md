---
title: "Useful matrix decompositions"
date: "2014-04-01"
questions:
- "How can matrix decompositions help with repeated solutions?"
- "What is the LU decomposition?"
- "What is the Cholesky decomposition?"
- "What is the QR decomposition?"
objectives:
- "Explain the main useful matrix decompositions"
katex: true
markup: "mmark"
---

Matrix factorisations play a key role in the solution of problems of the type $A x = b$. 
Often (e.g. ODE solvers), you have a fixed matrix $A$ that must be solved with many 
different $b$ vectors. A matrix factorisation is effectivly a pre-processing step that 
allows you to partition $A$ into multiple factors (e.g. $A = LU$ in the case of $LU$ 
decomposition), so that the actual solve is as quick as possible. Different 
decompositions have other uses besides solving $A x = b$, for example the $LU$, $QR$ and 
Cholesky decomposition (for a positive-definite matrix) can be used to to find the 
determinant of a large matrix. The Cholesky decomposition can be used to sample from a 
multivariate normal distribution, and is a very efficient technique to solve $A x = b$ 
for the specific case of a positive defintie matrix. The $QR$ decomposition can be used 
to solve a minimum least squares problem, to find the eigenvalues and eigenvectors of a 
matrix, and to calulcate the [Singular Value 
Decomposition](https://en.wikipedia.org/wiki/Singular_value_decomposition) (SVD), which 
is itself another very useful decomposition!

## $LU$ decomposition

The $LU$ decomposition is closely related to gaussian elimination. It takes the original 
equation to be solved $A x = b$ and splits it up into two separate equations involving a 
unit lower triangular matrix $L$, and the row echelon matrix $U$:

$$
\begin{aligned}
L y &= b \\
U x &= y
\end{aligned}
$$


where $A = LU$. The $L$ matrix is a *unit* lower triangular matrix and thus has ones on 
the diagonal, whereas $U$ is in row echelon form with pivot values in the leading 
coefficients of each row.

$$
A = \left( \begin{matrix}
1 & 0 & 0 & 0 \\ 
\ast & 1 & 0 & 0 \\
\ast & \ast & 1 & 0 \\
\ast & \ast & \ast & 1 \end{matrix}\right)
\left(\begin{matrix}
p_1 & \ast & \ast & \ast \\
0 & p_2 & \ast & \ast \\
0 & 0 & p_3 & \ast \\
0 & 0 & 0 & p_4
\end{matrix}\right)
$$

Thus, we have converted our original problem of solving $A x = b$ into two separate 
solves, first solving the equation $L y = b$, and then using the result $y$ to solve $U 
x = y$. 

$$
\begin{aligned}
\left(\begin{matrix}
1 & 0 & 0 & 0 \\ \ast & 1 & 0 & 0 \\
\ast & \ast & 1 & 0 \\
\ast & \ast & \ast & 1 
\end{matrix}\right)
\left(\begin{matrix}
y_1 \\ 
y_2 \\
y_3 \\
y_4 
\end{matrix}\right)
&= \left(\begin{matrix}
b_1 \\ 
b_2 \\
b_3 \\
b_4 
\end{matrix}\right)
\\
\left(\begin{matrix}
p_1 & \ast & \ast & \ast \\ 
0 & p_2 & \ast & \ast \\
0 & 0 & p_3 & \ast \\
0 & 0 & 0 & p_4 
\end{matrix}\right)
\left(\begin{matrix}
x_1 \\ 
x_2 \\
x_3 \\
x_4 
\end{matrix}\right)
&= \left(\begin{matrix}
y_1 \\ 
y_2 \\
y_3 \\
y_4 
\end{matrix}\right)
\end{aligned}
$$


However, each of those solves is very cheap to compute, in this case for the 4x4 matrix 
shown above the solution of $L y = b$ only needs 6 multiplication and 6 additions, 
whereas $U x = y$ requires 4 divisions, 6 multiplications and 6 additions, leading to a 
total of 28 arithmetic operations, much fewer in comparison with the 62 operations 
required to solve the original equation $A x = b$. In general, $LU$ decomposition for an 
$n \times n$ matrix takes about $2 n^3 / 3$ flops, or floating point operations, to 
compute.


## $LU$ factorisation without pivoting

A relativly simple $LU$ algorithm can be described if we assume that no pivoting is 
required during a gaussian elimination. In this case, the gaussian elimination process 
is a sequence of $p$ linear operations $E_1, E_2, ..., E_p$, with each operation $E_i$ 
being a row replacement that adds a multiple of one row to another below it (i.e. $E_i$ 
is lower triangular). The final matrix after applying the sequence of row reductions is 
$U$ in row echelon form, that is:

$$
E_p \cdots E_2 E_1 A = U
$$

Since we have $A = LU$, we can show that the sequence of operations $E_1, E_2, ..., E_p$ 
is also the sequence that reduces the matrix $L$ to an identity matrix:

$$
A = (E_p \cdots E_2 E_1)^{-1} U = LU,
$$

therefore, 

$$
L = (E_p \cdots E_2 E_1)^{-1},
$$

and,

$$
(E_p \cdots E_2 E_1) L = (E_p \cdots E_2 E_1) (E_p \cdots E_2 E_1)^{-1} = I
$$


This implies how we can build up the matrix $L$. We choose values for $L$ such that the 
series of row operations $E_1, E_2, ..., E_p$ convert the matrix $L$ to the identity 
matrix. Since each $E_i$ is lower triangular, we know that both $(E_p \cdots E_2 E_1)$ 
and $(E_p \cdots E_2 E_1)^{-1}$ are also lower triangular.

For example, consider the following matrix

$$
A = \left(\begin{matrix}
3 & 2 & 1 & -3 \\ 
-6 & -2 & 1 & 5 \\
3 & -4 & -7 & 2 \\
-9 & -6 & -1 & 15 
\end{matrix}\right)
$$

After three row reductions, $R_2 \mathrel{{+}{=}} 2 R_1$, $R_3 \mathrel{{+}{=}} -1 R_1$, 
and $R_3 \mathrel{{+}{=}} 3 R_1$, we have the following result:

$$
E_1 E_2 E_3 A = \left(\begin{matrix}
3 & 2 & 1 & -3 \\ 
0 & 2 & * & * \\
0 & -6 & * & * \\
0 & 0 & * & * 
\end{matrix}\right)
$$

To build the 1st column of $L$, we simply divide the 1st column of $A$ by the pivot 
value 3, giving

$$
L = \left(\begin{matrix}
1  & 0 & 0 & 0 \\ 
-2 & 1 & 0 & 0 \\
1  & * & 1 & 0 \\
-3 & * & * & 1 
\end{matrix}\right)
$$

For the next column we do the same, using the new pivot value $A_{2,2} = 2$ in row 2 to 
reduce $A_{3,2}$ and $A_{4,2}$ to zero, and then dividing the column vector under the 
pivot $(-6, 0)^T$ by the pivot value 2 to obtain the next column of $L$.

Repeating this process for all the columns in $A$, we obtain the final factorisation. 
You can verify for yourself that repeating the same row operations we did to form $U$ to 
the matrix $L$ reduces it to the identity matrix.

$$
L = \left(\begin{matrix}
1 & 0 & 0 & 0 \\ 
-2 & 1 & 0 & 0 \\
1 & -3 & 1 & 0 \\
-3 & 0 & 2 & 1 
\end{matrix}\right)
$$

$$
E_1 E_2 ... E_p A = U = \left(\begin{matrix}
3 & 2 & 1 & -3 \\ 
0 & 2 & 3 & -1 \\
0 & 0 & 1 & 2 \\
0 & 0 & 0 & 2 
\end{matrix}\right)
$$

## Pivoting

Of course, for any practial $LU$ factorisation we need to consider pivoting. Any matrix 
$A$ can be factorised into $PLU$, where $P$ is a permutation matrix, and $L$ and $U$ are 
defined as before. During the gaussian elimination steps we store an array of row 
indices $p_i$ indicating that row $i$ is interchanged with row $p_i$, and the resultant 
array of $p_i$ can be used to build the permutation matrix $P$ (It would be wasteful to 
store the entire martix $P$ so the array $p_i$ is stored instead). 

Thus, the LU algorithm proceeds as follows:

1. Begin with the left-most column $i=0$, find an appropriate pivot (e.g. maximum entry 
   in the colum) and designate this row as the pivot row. Interchange this row with row 
   $i$, and store the pivot row index as $p_i$. Use row replacements to create zeros 
   below the pivot. Create the corresponding column for $L$ by dividing by the pivot 
   value.
2. Continue along to the next column $i$, again choosing a pivot row $p_i$, 
   interchanging it with row $i$ and creating zeros below the pivot, creating the new 
   column in $L$, and making sure to record which pivot row has been chosen for each 
   column. Repeat this step for all the columns of the matrix.
3. Once the last column has been done, $U$ should be in row echlon form and $L$ should 
   be a unit lower triangular matrix. The array $p_i$ implicitly defines the permutation 
   matrix $P$

In practice, most library implementation store $L$ and $U$ in the same matrix since they 
are lower and upper triangular respectivly.

## Other Reading

- Linear algebra and its applications by David C. Lay. Chaper 2.5  
- Golub, G. H. & Van Loan, C. F. Matrix Computations, 3rd Ed. (Johns Hopkins University 
  Press, 1996). Chapter 3.2
- https://en.wikipedia.org/wiki/LU_decomposition

## Software

- 
  [`scipy.linalg.lu_factor`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.lu_factor.html).
- 
  [`scipy.linalg.lu_solve`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.lu_solve.html).
- 
  [`scipy.linalg.lu`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.lu.html).


## Problems

1. Take your gaussian elimination code that you wrote in the previous lesson and use it 
   to write an LU decomposition function that takes in a martix $A$, and returns $L$, 
   $U$ and the array $p_i$. You can check your answer using 
   [`scipy.linalg.lu_factor`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.lu_factor.html), 
   or by simply verifying that $PLU=A$
2. Write a unit test or set of unit tests using the Python `unittest` framework that you 
   can run to satisfy that your function is robust for a wide varietry of inputs $A$. 

## $LDL$ decomposition

It is often very benificial when solving linear systems to consider and take advantage 
of any special structure that the matrix $A$ might possesses. The $LDL$ decomposition is 
a varient on LU decomposition which is only applicable to a symmetric matrix $A$ (i.e. 
$A = A^T$). The advantage of using this decomposition is that it takes advantage of the 
redundent entries in the matrix to reduce the amount of computation to $n^3/3$, which is 
about a half that required for the $LU$ decomposition.

### Other reading


- Golub, G. H. & Van Loan, C. F. Matrix Computations, 3rd Ed. (Johns Hopkins University 
  Press, 1996). Chapter 4.1
- https://en.wikipedia.org/wiki/Cholesky_decomposition#LDL_decomposition_2

### Software

- 
[`scipy.linalg.ldl`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.ldl.html#scipy.linalg.ldl)

## Cholesky decomposition

*Symmetric positive definite* matrices are a very special type of matrix that often 
arise in practice. From a computational point of view, these class of matrix is very 
attractive because it is possible to decompose a symmetic positive definite matrix $A$ 
very efficiently into a single lower triangular matrix $G$ so that $A = GG^T$. 

A matrix $A$ is positive definite if $x^T A x > 0$  for any nonzero $x \in \mathbb{R}$. 
This statement by itself is not terribly intuative, so lets look at also look at an 
example of a $2 \times 2$ matrix

$$
A = \left(\begin{matrix}
a_{11} & a_{12} \\
a_{21} & a_{22}
\end{matrix}\right)
$$

If $A$ is symmetic positive definite (SPD) then

$$
\begin{aligned}
x &= (1, 0)^T \Rightarrow x^T A x = a_{11} > 0 \\
x &= (0, 1)^T \Rightarrow x^T A x = a_{22} > 0 \\
x &= (1, 1)^T \Rightarrow x^T A x = a_{11} + 2a_{12} + a_{22} > 0 \\
x &= (1,-1)^T \Rightarrow x^T A x = a_{11} - 2a_{12} + a_{22} > 0 \\
\end{aligned}
$$

The first two equations show that the diagonal entries of $A$ must be positive, and 
combining the last two equations imply $|a_{12}| \le (a_{11} + a_{22}) / 2$, that is 
that the matrix has much of its "mass" on the diagonal (note: this is *not* the same as 
the matrix being diagonally dominant, where $|a_{ii}| > \sum_{i=1...n,j \ne i} 
|a_{ij}|$). These two observations for our $2 \times 2$ matrix also apply for a general 
$n \times n$ SPD matrix. One of the very nice consequences of this "weighty" diagonal 
for SPD matrices is that it precludes the need for pivoting.

It can be shown that if $A$ is a SPD matrix, then the $LDL^T$ decomposition exists and 
that $D = \text{diag}(d_1, ..., d_n)$ has positive diagonal entries. Therefore, it is 
straightforward to see that $LDL^T$ = $GG^T$, where $G = L \text{diag}(\sqrt{d_1}, ..., 
\sqrt{d_n})$. The decomposition $A = GG^T$ is known as the cholesky decomposition and 
can be efficiently constructed in $n^3 / 3$ flops. There are a number of algorithms to 
construct this decomposition, and both the [wikipedia 
entry](https://en.wikipedia.org/wiki/Cholesky_decomposition) and Chapter 4.2 of the 
Matrix Computations textbook by Golub and Van Loan gives a number of different varients.

Note that a $LDL$ decomposition can also be used to calculate a cholesky decomposition, 
and this could be more efficient approach since (a) the SPD structure means that we can 
neglect pivoting in the $LDL$ decomposition, and (b) the $LDL$ decomposition does not 
requiring taking the square root of the diagonal elements. 

### Other Reading

- Golub, G. H. & Van Loan, C. F. Matrix Computations, 3rd Ed. (Johns Hopkins University 
  Press, 1996). Chapter 4.2
- https://en.wikipedia.org/wiki/Cholesky_decomposition

### Software

- 
  [`scipy.linalg.cholesky`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cholesky.html)
- 
[`scipy.linalg.cho_factor`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cho_factor.html)
- 
[`scipy.linalg.cho_solve`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cho_solve.html)
- 
[`scipy.linalg.cholesky_banded`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cholesky_banded.html)
- 
[`scipy.linalg.cho_solve_banded`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cho_solve_banded.html)

## Problems

Imagine that we wanted to sample an array of values $x_i$, for $i = 1...n$, where each 
value is sampled from an independent Normal distribution with standard deviation 
$\sigma$

 $$x_i \sim \mathcal{N}(0, \sigma)$$

 This could be achieved, for example, by sampling from a Normal distribution with unit 
 standard deviation, a function that typically exists in any computer language, then 
 multiplying by $\sigma$

 $$x_i = \sigma \eta$$

 where $\eta \sim \mathcal{N}(0, 1)$

 Now imagine that instead of an independene Normal distribution you wish to sample 
 $\mathbf{x} = [x_1, x_2, ..., x_n]$ from a multivariate Normal distribution with some 
 covariance matrix $\Sigma$

 $$\mathbf{x} \sim \mathcal{N}(\mathbf{0}, \Sigma)$$

 We can achive this in practice by using the Cholesky decomposition. A covariance 
 matrix is a symmetic positive semidefinite matrix (i.e. $x^T \Sigma x \ge 0$}, and 
 therefore can be decomposed into  $\Sigma = LL^T$. We can then draw a sample from 
 $\mathcal{N}(\mathbf{0}, \Sigma)$ by scaling an independently generated random vector 
 by $L$

 $$\mathbf{x} = L \mathbf{\eta}$$

 where each element of the vector $\eta$ is $\eta_i \sim \mathcal{N}(0, 1)$.

 Write Python code to randomly sample an n-dimensional vector $x$ from 
 
 1. an indepenindependent Normal distribution with variance $\sigma^2$

 2. a multivariate normal distribution using a covariance matrix $\Sigma_{ij} = \exp[(i- 
    j)^2/ \sigma^2]$

 3. a multivariate normal distribution with $\Sigma = \sigma^2 I$. Show that this 
    algorithm reduces to that used for (1).

## QR decomposition

## The least-squares problem

One of the most important application of the $QR$ decomposition is the least squares 
solution of a set of overdetermined equations. That is a set of $m$ linear equations 
with $n$ unknowns, with $m \ge n$. The least squares problem to be solved is the 
mimimisation of $||A x - b ||_2$, where $|| x ||_2 = \sqrt{x_1^2 + x_2^2 + ... + x_m^2}$ 
is the standard 2-norm, and where $A \in \mathbb{R}^{m \times n}$ with $m \ge n$ and $b 
\in \mathbb{R}^m$. In this case, the problem $Ax = b$ will often have no solution, and 
thus it is nessessary to consider $Ax$ and $b$ as *approximatelly* equal, and to 
minimise the distance between them by minimising the loss function $||A x - b||_2$.

To solve this least squares problem, we need to consider the subspace of all vectors in 
$\mathbb{R}^{m}$ that are formed from linear combinations of the columns of $A$. This is 
known as the column space of the $A$, and is denoted as Col $A$. Given that *any* linear 
combination of the columns of $A$ will lie in this space, we can say that $Ax$ will also 
lie in Col $A$ for any $x$.

Now consider a projection of $b$ into the column space of $A$ to give a new vector 
$\hat{b}$ (i.e. $\hat{b}$ is the closest point in Col $A$ to $b$), see the diagram 
below. Because $\hat{b}$ is in the column space of $A$, we know that there is another 
vector $\hat{x}$ that also lies in Col $A$ and satisfies

$$
A \hat{x} = \hat{b}
$$

Since $\hat{b}$ is the closest point to $b$ in the column space of $A$, we can therefore 
say that $\hat{x}$ is the least-squares solution.


![least squares problem](/figs/linear-least-squares.svg)


We can show that the vector $b - \hat{b} = b - A \hat{x}$ is orthogonal to Col $A$ and 
therefore also orthogonal to each column in $A$, so we have $a_j^T (b - A \hat{x})$ for 
each column $a_j$ of $A$. Putting these $m$ equations together we can write

$$
A^T (b - A \hat{x}) = 0
$$

or rearranged slightly, we can find the least-sqaures solution $\hat{x}$ via the 
solution of the equation

$$
A^T A \hat{x} = A^T b
$$

The $QR$ decomposition divides $A = QR$ into an orthogonal matrix $Q$, and an upper 
triangular matrix $R$. Most importantly for the least-squares problem, the matrix $Q$ is 
also an orthonormal basis for Col $A$ and therefore $\hat{b} = Q Q^T b$.

Given this decomposition, it can be shown that the least squares solution of $A x = b$ 
is given by

$$
\hat{x} = R^{-1} Q^T b
$$

To prove this, we let $\hat{x} = R^{-1} Q^T b$ and make the following substitutions

$$
A\hat{x} =  QR \hat{x} = QRR^{-1}Q^T b = Q Q^T b = \hat{b}
$$

Therefore $A\hat{x} = \hat{b}$, which proves that $\hat{x}$ is the least-squares 
solution for $A x = b$

Finally, we note that the inverse $R^{-1}$ should not be calculated directly, but 
instead $\hat{x}$ should be found by solving

$$
R x = Q^T b
$$

### Constructing the QR decomposition

$QR$ decomposisions are normally computed via Householder reflections, Givens rotations 
or the Gram-Schmidt process. For a brief summary of the first two methods, it is useful 
to consider a simple $2 \times 2$ reflection or rotation of a 2d vector. For example, 
the matrix

$$
Q = \left(\begin{matrix}
\cos(\theta) & \sin(\theta) \\
-\sin(\theta) & \cos(\theta)
\end{matrix}\right)
$$

is a *rotation* matrix that when applied to a vector $x$ will result in $y = Qx$, where 
$y$ is rotated counterclockwise through the angle $\theta$. $Q$ is also *orthogonal* 
since $QQ^T = I$.

Similarly, a $2 \times 2$ *reflection* matrix can be constructed as 

$$
Q = \left(\begin{matrix}
\cos(\theta) & \sin(\theta) \\
\sin(\theta) & -\cos(\theta)
\end{matrix}\right)
$$

which when applied to a vector $x$ will result in $y = Qx$, where $y$ is reflected 
across the line defined by $\text{span}((\cos(\theta), \sin(\theta))^T)$.

Rotations and reflections are often useful because they can be selected in order to 
introduce zeros to the vector they are applied to. Given an $m \times n$ matrix $A$, a 
series of $n$ *Householder reflections* can be applied to reduce $A$ to an upper 
triangular matrix $R$

$$
H_n ... H_2 H_1 A = R
$$

By setting $Q = H_1 H_2 ... H_n$, we can show that $A = QR$, and that $Q$ is an 
orthogonal matrix which is also an orthonormal basis for the column space of $A$.

Similarly, a *Givens rotation* can be used to zero a single component of $A$, so that a 
a series of rotations can be used to contruct the upper triangular matrix $R$

$$
G_j ... G_2 G_1 A = R
$$

so that $Q = G_1 G_2 ... G_j$, and $A = QR$. For both the Householder and Givens 
methods, it is often useful to not construct the full matrix $Q$ but to keep $Q$ 
factored as a implicit product of either $H_1 H_2 ... H_n$ or $G_1 G_2 ... G_j$. Fast 
algorithms exist to calculate the produce of these factored forms to another vector.

The final method to contruct a $QR$ decomposition is using the Gram-Schmidt process, 
which is a process for contructing an orthogonal or orthonormal basis for a given 
subspace defined by the span of the set of vectors $x_1, x_2, ..., x_n$. If these $n$ 
vectors are the columns of the $m \times n$ matrix $A$, then the Gram-Schmidt process 
can be used to directly contruct the orthonormal basis of the column space of $A$ given 
by $Q$, and that $A = QR$ where $R$ is an upper triangular matrix. The matrix $R$ can be 
calculated using $R = Q^T A$. Note that the classical Gram-Schmidt exhibits poor 
numerical qualities, therefore a modified version of the algorithm exists, which is 
described in the Golub and Van Loan Matrix Computations textbook listed below.

In terms of computational work, the Householder method takes $2n^2(m-n/3)$ flops to 
compute $Q$ in factored form, and another $2n^2(m-n/3)$ to get the full matrix $Q$, 
whereas the Gram-Schmidt method is more efficient at $2mn^2$ flops. However, Householder 
is normally prefered in practice as even with the modified algorithm the numerical 
properies of the Gram-Schmidt are still poor in comparison with both Householder and 
Givens (i.e. the final orthogonality of $Q$ is not ideal), so is only useful when the 
columns of $A$ are already fairly independent. Using Givens rotations the matrix $R$ can 
be found in $2n^2(m-n/3)$, or the factorised form of the $QR$ decomposition can be found 
in the same amount of time. The full matrix $Q$ is not normally calculated via Givens 
rotations. Using Givens rotations is most useful when there are only few non-zeros in 
$A$, and is more easily parallised than Householder.

### Other Reading

The discussion in this section relied on concepts such as orthogonal and orthonormal 
vector pairs, vector spaces and subspaces and basis vectors. It is well worth 
investigating these topics further in:

- Linear algebra and its applications by David C. Lay. Chapers 4 & 6.

Additional reading on the $QR$ decomposition can be found at:

- Linear algebra and its applications by David C. Lay. Chaper 6.4  
- Golub, G. H. & Van Loan, C. F. Matrix Computations, 3rd Ed. (Johns Hopkins University 
  Press, 1996). Chapter 5.2
- https://en.wikipedia.org/wiki/QR_decomposition

### Software

- 
  [`scipy.linalg.qr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.qr.html)
- 
  [`scipy.linalg.qr_multiply`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.qr_multiply.html)
- 
  [`scipy.linalg.qr_update`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.qr_update.html)
- 
  [`scipy.linalg.qr_delete`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.qr_delete.html)
- 
  [`scipy.linalg.qr_insert`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.qr_insert.html)


## Problems

For this exercises we will be using some data on Oxford's weather which is hosted by 
[Saad Jbabdi](https://users.fmrib.ox.ac.uk/~saad/) from the Wellcome Centre for 
Integrative NeuroImaging (FMRIB), which can be obtained 
[here](http://www.fmrib.ox.ac.uk/~saad/ONBI/OxfordWeather.txt).

We wish to fit a quadratic model of the form $y = a x^2 + b x + c$ to the hours of 
sunlight observed in Oxford (7th column in `OxfordWeather.txt`) versus the month (2nd 
column). The dataset in question has $m > 3$ data points, so our model gives us a set of 
$m$ equations for 3 unknowns $a$, $b$, and $c$ that are overdetermined, that is, for 
each data point $(y_i, x_i)$ for $i=1..m$ we have:

$$
y_i = a x_i^2 + b x_i + c
$$

Use a $QR$ decomposition to find the least-squares solution to these equations, and 
therefore fit the model to the data. Plot the model and the data side by side to 
qualitativly evaluate the fit.

