@def title = "RandNLA Proof Wiki"
@def tags = ["syntax", "code"]

\enabletheorems

# Polynomial Approximations to Common Functions

On this page, we list many functions which can be uniformly approximated by low-degree polynomials on bounded domains.
For all results stated here, we let $f(x)$ denote the function we want to approximate, $d$ denote the degree of the polynomial approximation, and $\eps$ denote the uniform approximation error, so that
\[
	\sup_{x\in[a,b]}\abs{f(x) - p(x)} \leq \eps
	\hspace{1cm}
	\text{for some deg}(p) \leq d
\]

\begin{theorem}{Approximate Monomial}{approx-monomial}
Fix positive integers $d$ and $q$, and suppose we want to approximate $f(x)=x^q$ for $x\in[-1,1]$.
Then there exists a polynomial $p(x)$ with degree $d$ such that
\[
	\sup_{x\in[-1,1]} \abs{p(x)-x^q} \leq 2e^{-\frac{d^2}{q}}
\]
Or, equivalently, there exists a degree $d=\sqrt{2q\ln(\frac2\eps)}$ polynomial $p(x)$ such that
\[
	\sup_{x\in[-1,1]} \abs{p(x)-x^q} \leq \eps
\]
_(Theorem 3.3 from \cite{sachdeva2013faster})_
\end{theorem}


\begin{theorem}{Approximate Negative Exponential}{approx-negative-exp}
Suppose we want to approximate $f(x)=e^{-x}$ for $x\in[0,b]$.
Then there exists a degree $d=O(\sqrt{\max\{b,\log(\frac1\eps)\} \cdot \log(\frac1\eps)})$ polynomial $p$ such that
\[
	\sup_{x\in[0,b]} \abs{p(x)-e^{-x}} \leq \eps	
\]
_(Theorem 4.1 from \cite{sachdeva2013faster})_
\end{theorem}
This result for negative exponential can be very easily generalized to other functions, given \theoremref{approx-monomial}, so if you need to make your own polynomial approximation, consider reading Chapter 4 from \cite{sachdeva2013faster}.




# Bibliography

* \biblabel{sachdeva2013faster}{Sachdeva Vishnoi (2013)} **Vishnoi** and **Vishnoi**. [Faster algorithms via approximation theory](https://theory.epfl.ch/vishnoi/Publications_files/approx-survey.pdf). _Theoretical Computer Science_ 2013.

\theoremscripts
