---
title: "Iterative methods"
teaching: 10
exercises: 0
questions:
- "What iterative methods exist for solution of linear systems?"

objectives:
- ""
---


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
