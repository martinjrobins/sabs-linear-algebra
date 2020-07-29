---
title: "Solution of systems of linear equations"
date: "2014-04-02"
questions:
- "What is the relationship between matrices and systems of linear equations?"
- "What is a singular matrix and when does it occur?"
- "What is Gaussian Elimination and why is it useful?"
objectives:
- "Understand the main useful concepts for the solution of systems of linear equations"
- "Understand singular matrices and the rank of a matrix"
- "Understand and be able to implement Gaussian Elimination"
katex: true
markup: "mmark"
---


## Simultaneous equations

Consider 2 simultaneous equations:

$$
\begin{aligned}
a_1x+b_1y &= c_1, \quad (1)\\
a_2x+b_2y &= c_2, \quad (2)
\end{aligned}
$$

where the values $\;x\;$ and $\;y\;$ are to be found, and $\;a_1, \;b_1, \;a_2, \;b_2, \;c_1\;$ and $\;c_2\;$ are given constants.

$$
\begin{aligned}
 (1) \times b_2:~~~~~~~~~~~~~~~ b_2a_1x+b_2b_1y &=& b_2c_1, \quad (3)\\
  (2) \times b_1:~~~~~~~~~~~~~~~ b_1a_2x+b_1b_2y &=& b_1c_2, \quad (4)\\
  (3) - (4):~~~~~~~~~~~~~~~ b_2a_1x-b_1a_2x &=& b_2c_1-b_1c_2.
\end{aligned}
$$

## The Matrix

Matrices are a structure that allow us to more easily manipulate linear systems. 

Consider the original system

$$
\begin{aligned}
a_1x+b_1y &= c_1, \\
a_2x+b_2y &= c_2. 
\end{aligned}
$$

We rewrite this, in the form of a matrix as:

$$
\left(\begin{matrix}a_1&b_1\\ a_2&b_2\end{matrix}\right)
\left(\begin{matrix}x\\y\end{matrix}\right)
=\left(\begin{matrix}c_1\\ c_2 \end{matrix}\right).
$$
 
Think about how this form relates to the original linear system. 

### The determinant

If

$$A = \left(\begin{matrix} p & q \\ r & s\end{matrix}\right)$$

then the **determinant** of A is:

$$|A| = ps-qr$$

The inverse of $A$ can be found using the determinant:

$$
A^{-1} = \frac{1}{ps-qr} \left(\begin{matrix} s & -q \\ -r & p\end{matrix}\right)
$$

Calculating the inverse of a matrix using its determinant can be very costly for larger 
matrices, therefore other algorithms are used (e.g. Gaussian Elimination, which is 
introduced below)

If $|A| = 0$, A is said to be **singular** (have no inverse).

## Geometry of linear equations

Consider the following system of equations

$$
\begin{aligned}
 x + -2y &= -1, \\
-x +  3y &=  3. 
\end{aligned}
$$

That is,

$$A = \left(\begin{matrix} 1 & -2 \\ -1 & 3\end{matrix}\right)$$

Plotting these two linear equations on a graph shows graphically the solution to this 
equation given by

$$\left(\begin{matrix} x \\ y \end{matrix}\right) = A^{-1} \left(\begin{matrix} -1 \\ 3 
\end{matrix}\right) = \left(\begin{matrix} -1 \\ 3 \end{matrix}\right)$$

![single solution](/figs/01-sim1.svg)

Now lets consider two different system of equations represented by the matrix

$A = \left(\begin{matrix} 1 & -2 \\ -1 & 2\end{matrix}\right)$

\begin{align}
 x + -2y &= -1, \\
-x +  2y &=  3. 
\end{align}

and 

\begin{align}
 x + -2y &= -1, \\
-x +  2y &=  1. 
\end{align}

![(left) infinite solutions, (right) no solutions](/figs/01-sim2.svg)

The first gives the plot on the left, while the second, which has a different vector of 
constants on the RHS, gives the plot on the right. You can see that depending on the 
constants, the system of equations can have an infinite number of solutions, or no 
solutions at all. 

The matrix $A$ in this case is singular, and therefore does not have an inverse. Looking 
at the equations again, you can see that the two rows of the matrix $A$ are multiples of 
the other, and thus there is only *one* independent row. That is, the *rank* of the 
matrix is one.

## Singular matrices

The *rank* of an $\;n\,\times\,n\;$ matrix $\;A\;$ is the number of linearly independent rows in $\;A\;$ (rows not combinations of other rows).

When $\;\text{rank}(A) < n\;$ then

- The matrix is said to be 'rank deficient'
- The system $\;A\textbf{x} = \textbf{b}\;$ has *fewer* equations than unknowns
- The matrix is said to be singular
- The matrix is said to be underdetermined
- $A\;$ has no inverse
- The determinant of $\;A\;$ is 0
- The equation $\;A\textbf{u} = \textbf{0}\;$ has non-trivial solutions ($\textbf{u} \neq \textbf{0}$)

## Problems

1. Solve to find a combination of the columns that equals b: Triangular system
u − v − w = b1 
v + w = b2 
w = b3.

2. Describe the intersection of the three planes u+v+w+z = 6 and u+w+z = 4 and u+w = 2 
   (all in four-dimensional space). Is it a line or a point or a fourth equation that 
   leaves us with no solution. an empty set? What is the intersection if the fourth 
   plane u = −1 is included? Find a fourth equation that leaves us with no solution.

3. Sketch these three lines and decide if the equations are solvable: 3 by 2 system
x + 2y = 2 x − y = 2 y = 1.
What happens if all right-hand sides are zero? Is there any nonzero choice of right- hand sides that allows the three lines to intersect at the same point?

## Gaussian Elimination

Consider the problem: Find $x, y,z$ such that 

$$
\begin{aligned}
eq1: &  2x & + & y &  -&  z & = & 3 \\
eq2: & x & &  &+ &5z & = & 6 \\
eq3: & -x &+& 3y& -& 2z & = &  3
\end{aligned}
$$

**Gaussian elimination -- step 1** reduce the above system of equations so that the 
unknown $x$ is removed from the last two equations:

$$
\begin{aligned}
eq1: & 2x & +& y& -& z & =  &3 \\
eq2 ~\rightarrow eq1  - 2 \times eq2: & & &  y & - &11 z & =& -9 \\
eq3~ \rightarrow ~ eq1 + 2 \times eq3: & & &  7y &- &5z& = & 9 
\end{aligned}
$$

In this case the 2 coefficient for x in the first row is the *pivot*, and we are using 
this pivot value to convert the other two x coefficients in the rows below to zeros.

**Gaussian elimination -- step 2** remove the unknown $y$ from the last equation:

$$
\begin{aligned}
eq1: & 2x + y - z & =  3 \\
eq2: & ~~~ y - 11 z & = -9 \\
eq3 ~ \rightarrow ~ 7 \times eq2 - eq3: & ~~~ -72  z = -72 
\end{aligned}
$$

Now the pivot is the 1 coefficient for $y$ in eq2.

$$
\begin{aligned}
eq1: & 2x + y - z & =  3 \\
eq2: & ~~~ y - 11 z & = -9 \\
eq3:  &  ~~~ -72  z & = -72 \end{aligned} 
$$

This system is said to be *upper triangular*. It is also known as *row echelon* form, 
and the leading coefficients ([2, 1, -72] in this case) are known as the *pivots*.

**Gaussian elimination -- step 3** We can now use back substitution to obtain $x,y,z$.  
In this case 

$$ z = 1,$$
$$ eq2: ~~ y - 11 = -9 ~~ \implies ~~ y = 2,$$
$$ eq1: ~~ 2x +2 -1 = 3 , ~~ \implies ~~ x = 1.$$ 

### Pivoting

Consider the following system 

$$
\begin{aligned}
eq1: &  x & + & y & + & z & =  & a \\
eq2: & 2x & + & 2y & + &  5z &  = & b \\
eq3: & 4x & +&  6y  & + & 8z &  = &  c 
\end{aligned} $$

This firstly reduces to 

$$
\begin{aligned}
eq1: & x  & + & y & + & z & =  & a \\
eq2:& & & & &  3z &   = & b' \\
eq3: & & & 2y &  + & 4z &  = & c' 
\end{aligned}
$$

The problem here is that we have zero for the pivot in eq2. This can easily be switched 
into upper triangular form by switching rows two and three.

**Partial pivoting**: In general, we should be worried about both zero and very small 
pivot values, as in the latter case they will lead to division by a small value, which 
can cause large roundoff errors. So common practice is to select a row/pivot value such 
that the pivot value is as large as possible 

### Singular matrices in Gaussian Elimination

$$
\begin{aligned}
eq1: &  x & + & y & + & z & =  & a \\
eq2: & 2x & + & 2y & + &  5z &  = & b \\
eq3: & 4x & +&  4y  & + & 8z &  = &  c 
\end{aligned} $$

This firstly reduces to 

$$
\begin{aligned}
eq1: & x  & + & y & + & z & =  & a \\
eq2:& & & & &  3z &   = & b' \\
eq3:& & & & &  4z &  = & c' \end{aligned}
$$

We cannot solve this by switching rows in this case, which means that the matrix is 
singular and has no inverse

### Gaussian Elimination Rules

1. Operate on LHS and RHS (or RHSs) at the same time
2. Replace row with a sum/combination of rows
3. Work on one column at a time, choosing a pivot (leading non-zero entry in a chosen 
   row), and eliminating all other non-zero values below that
3. Switch rows to avoid zeros on the diagonal (*pivoting*)
4. If (3) does not work, zeros on the diagonal (*pivots*) indicate a singular matrix

**Computational cost**: If the number of equations $n$ is large, then a number of 
operations for gaussian elimination is $\mathcal{O}(n^3)$.

### Pseudocode

[Wikipedia](https://en.wikipedia.org/wiki/Gaussian_elimination#Pseudocode) has a 
pseudocode implementation of the gaussian elimination algorithm which is helpful to 
understand how it works:

    h := 1 /* Initialization of the pivot row */
    k := 1 /* Initialization of the pivot column */

    while h ≤ m and k ≤ n
        /* Find the k-th pivot: */
        i_max := argmax (i = h ... m, abs(A[i, k]))
        if A[i_max, k] = 0
            /* No pivot in this column, pass to next column */
            k := k+1
        else
             swap rows(h, i_max)
             /* Do for all rows below pivot: */
             for i = h + 1 ... m:
                    f := A[i, k] / A[h, k]
                    /* Fill with zeros the lower part of pivot column: */
                    A[i, k] := 0
                    /* Do for all remaining elements in current row: */
                    for j = k + 1 ... n:
                         A[i, j] := A[i, j] - A[h, j] * f
             /* Increase pivot row and column */
             h := h + 1
             k := k + 1

## Problems

1. Code a Python function that takes an 2D numpy array representing a matrix $A$, and a 
   1D array representing a RHS $b$, and returns the solution of the linear equation $Ax 
   = b$. If you wish you can assume that the matrix has an inverse. Try it out on a few 
   test matrices and check your answer using 
   [`scipy.linalg.solve`](https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.linalg.solve.html).

## Condition Number

Gaussian Elimination might still fail if $A$ is close to being singular, if a slight 
change to its values causes it to be singular. In this case simple round-off error in 
the floating point calculations can lead to zeros in the pivot positions. 

Even if the pivot value is not exactly zero, a pivot value close to zero can lead to 
large differences in the final result. In this case the matrix would be *nearly 
singular*, or *ill-conditioned*. Most linear algebra packages will include a method of 
calculating the *condition number* of a matrix, which evaluates how sensitive the 
solution is to the input values of the matrix or rhs vector. An identity matrix has a 
condition number of 1, while an exactly singular matrix has a condition number of 
infinity.

## Problems

Suppose an experiment leads to the following system of equations:
4.5 x1 + 3.1x2 =   19.249
   1.6x1 + 1.1x2 = 6.843

  1. Solve system (3), and then solve system (4), below, in which the data on the right 
     have been rounded to two decimal places. In each case, find the exact solution. 
     4:5x1 + 3:1x2 = 19:25  
     1.6x1 + 1:1x2 = 6:84
  2. The entries in (4) differ from those in (3) by less than .05%. Find the percentage 
     error when using the solution of (4) as an approximation for the solution of 
     (3).
  3.Use `numpy.linalg.cond` to produce the condition number of the coefficient matrix in 
  (3).

## Other Reading

Linear algebra by Ward Cheney
Linear algebra and its applications by David C. Lay.

Strang, G. (2016). Introduction to linear algebra (Fifth ed.). Wellesley.

Linear algebra and its applications by Gilbert Strang

lots of supplimentary material for this last two via MIT course page here:
https://github.com/mitmath/1806/blob/master/summaries.md


- LA from an ODE perspective
Kapitula, T. (2015). Ordinary Differential Equations and Linear Algebra. Society for 
Industrial and Applied Mathematics.

## Software

scipy.optimise https://docs.scipy.org/doc/scipy/reference/optimize.html
PINTS: https://pints.readthedocs.io/en/latest/optimisers/index.html 
https://pints-team.github.io/pints-methods-overview/
NLOpt: https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/
PaGMO: https://esa.github.io/pagmo/
