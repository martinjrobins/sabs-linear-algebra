---
title: "Sparse matrices"
teaching: 10
exercises: 0
questions:
- "How can matrix decompositions help with repeated solutions?"
- "What is the LU decomposition?"
- "What is the QR decomposition?"
- "What is the Cholesky decomposition?"

objectives:
- "Explain what a sparse matrix is"
- "Demonstrate how sparse matices can lead to speed-up for realistic matrices"
---

## Problems

Do these but using scipy.sparse instead of matlab:

item {\it Sparse matrices}

For many of the matrices you will want to manipulate in {\sc Matlab} most of the entries will be zero. This could be protein interaction data or could be result 
of numerical approximation of partial differential equations. These types of matrices are known as sparse. {\sc Matlab} has sparse versions of most of the inbuilt functions, which operate much faster than their dense counterparts. 

\begin{enumerate}
\item Enter the following commands into {\sc Matlab} and note what they do and what you see. You may need to use the \mcode{doc} and \mcode{help} commands.
 
\begin{lstlisting}
	A=diag(diag(2*ones(100)))+diag(ones(99,1),1)+diag(ones(99,1),-1);
	nnz(A)
	spy(A)
	Ainv=inv(A);
	spy(Ainv)
	B=sparse(A);
	spy(B)
	spy(A-B)
\end{lstlisting}
Use the command \mcode{whos} to see which matrix, \mcode{A} or \mcode{B}, uses the most memory. Why does \mcode{Ainv} have lots of non zero entries?

The following commands are useful in building sparse matrices enter them into {\sc Matlab} to see what they do.
\begin{lstlisting}
	speye(5)
	C=sprand(10,10,0.1)
	[i,j,v]=find(C)
	[i,j,v]=find(C~=0)
	[i,j]=find(C>1)
	
    spdiags(ones(10,1),[1],10,10)
	spdiags(ones(10,1)*[-1,2,-1],[-1,0,1],10,10)
	clear
\end{lstlisting}

\item Using the commands \mcode{tic} and \mcode{toc} we can time how long it takes to execute a command. Enter the following commands
into {\sc Matlab}.
\begin{lstlisting}
  A=2*eye(500)-diag(ones(499,1),1) ...
      -diag(ones(499,1),-1);
  tic, A*A; toc

  B=spdiags(ones(500,1)*[-1,2,-1],[-1,0,1],500,500);
  tic, B*B; toc
\end{lstlisting}
Which takes longer?
\item Use the  \mcode{tic} and \mcode{toc} commands to compare the efficiency of the ``\mcode{\\}'' (backslash) operator, on dense and sparse matrices using \mcode{A} and \mcode{B} from above and the following vector.
\begin{lstlisting}
  b=ones(500,1);
\end{lstlisting}
\end{enumerate}


The remaining exercises involve the manipulation and solution of the linear system 
resulting from the finite difference solution to Poisson's equation.
Let $A$ be a sparse symmetric positive definite matrix of dimension $(N-1)^2
\times (N-1)^2$ entered in \textsc{Matlab} (for a given $N$) by the function
\mcode{buildA.m} as follows:
\vspace*{-1.5mm}
\begin{lstlisting}
  function A=buildA(N)
  dx=1/N;
  nvar=(N-1)^2;
  e1=ones(nvar,1);
  e2=e1; e2(  1:N-1:nvar)=0;
  e3=e1; e3(N-1:N-1:nvar)=0;
  A=spdiags([-e1 4*e1 -e1],-(N-1):N-1:N-1,nvar,nvar) ...
   +spdiags([-e3      -e2],    -1: 2 :  1,nvar,nvar);
  A=A/dx^2;
\end{lstlisting}
\vspace*{-1.5mm}
and let $\mathbf{f}_1$ and $\mathbf{f}_2$ be the vectors defined in
\mcode{buildf1.m} and \mcode{buildf2.m}
\vspace*{-1.5mm}
\begin{lstlisting}
  function f1=buildf1(N)
  x=0:1/N:1; y=x;
  f=sin(pi*x)'*sin(pi*y);
  f1=reshape(f(2:N,2:N),(N-1)^2,1);
\end{lstlisting}
\vspace*{-2.5mm}
\begin{lstlisting}
  function f2=buildf2(N)
  x=0:1/N:1; y=x;
  f=max(x,1-x)'*max(y,1-y);
  f2=reshape(f(2:N,2:N),(N-1)^2,1);
\end{lstlisting}
\vspace*{-1.5mm} 
The files \mcode{buildA.m}, \mcode{buildf1.m} and \mcode{buildf2.m} are all available on web learn.
We will consider manipulation of the matrix $A$ and solution of the linear
systems $A\mathbf{U}_i=\mathbf{f}_i$. The solution to this linear system
corresponds to a finite difference solution to Poisson's equation $-\nabla^2 u
= f$ on the unit square with zero Dirichlet boundary conditions where $f$ is
either $\sin(\pi x) \sin (\pi y)$ or $\max(x,1-x) \max(y,1-y)$. As
you'll discover later in the year, PDEs of this type occur (usually
with some additional reaction and or convection terms) very frequently
in mathematical modelling of physiological processes, and even in image
analysis. 

\begin{enumerate}
\setcounter{enumi}{7}


\item {\it Vector and matrix manipulation}

\begin{enumerate}
\item Set $N=4$ and produce a spy plot of the matrix $A$. How many nonzero
  diagonals are there? How many non-zero entries are there?
\item What is the determinant of $A$ when $N=4$?
\item When $N=4$ which entries of $A^{-1}$ are greater than $0.02$?
\item For $N=4$ check that $A$ is symmetric by producing a spy plot of
  $A-A^T$: if $A$ is symmetric there should be no non-zero entries. Check
  that $A$ is positive definite by looking at the eigenvalues of $A$: they
  should all be strictly positive.
\end{enumerate}


