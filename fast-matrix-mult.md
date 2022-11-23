@def title = "RandNLA Proof Wiki"
@def tags = ["syntax", "code"]

\enabletheorems
\newcounter{NumAlgorithms}

# Fast Matrix Multiplication

On this page we show one of several possible fast matrix multiplication guarantees, where we approximate $\mA^\intercal\mB$ as $(\mS\mA)^\intercal(\mS\mB)$, where is a randomized sampling matrix.
Notably, these guarantees generally have a polynomial dependence on the success probability, though they can be boosted to have a $\log(\frac1\delta)$ success probability by repeating the algorithm.

We consider the following algorithm:
\begin{algorithm}{Fast Matrix Multiplication}{fast-mmult}
**input**: Matrices $\mA\in\bbR^{n \times d}$, $\mB\in\bbR^{n \times m}$. Number $k$ of subsamples. Probabilities $p_1,\ldots,p_n$

**output**: Sketched matrix $\mC\in\bbR^{d \times m}$

1. Sample indices $s_1,\ldots,s_k\in[n]$ iid wrt $p_1,\ldots,p_n$
1. Build the sample-and-rescale matrix $\mS\in\bbR^{k \times n}$:
    
    Row $t$ of $\mS$ has form $\begin{bmatrix}0&0&\cdots&0&\frac{1}{\sqrt{k p_{s_t}}}&0&\cdots&0\end{bmatrix}$, where index $s_t$ is the nonzero entry.

1. Return $\mC = (\mS\mA)^\intercal(\mS\mB)$
\end{algorithm}


Since $(\mS\mA)^\intercal \in \bbR^{d \times k}$ and $\mS\mB\in\bbR^{k \times m}$, we can compute $(\mS\mA)^\intercal(\mS\mB)$ in $O(kdm)$ time instead of $O(ndm)$ time, just using naive matrix multiplication.
We show the following:

\begin{theorem}{Fast Matrix Multiplication}{fast-mmult}
Fix $\eps>0$ and $\delta\in(0,1)$.
Let $\mC$ be the resulting of fast matrix multiplication with $k \geq \frac{1}{\eps^2\delta}$ and $p_\ell = \frac{\normof{\va_\ell}_2^2}{\normof{\mA}_F^2}$, where $\va_\ell$ is the $\ell^{th}$ row of $\mA$.
Then with probability $1-\delta$,
\[
	\normof{\mC-\mA^\intercal\mB}_F \leq \eps\normof{\mA}_F \normof{\mB}_F
\]
\end{theorem}

Notably, we are not hiding any constants when we say $k \geq \frac1{\eps^2 \delta}$.
We prove the results in two steps.
First, we show a result for arbitrary sampling probabilities, then we prove \theoremref{fast-mmult}.
Also, there's a lot of indexing in this analysis, so to be clean we consistently denote $t \in [k]$, $\ell \in [n]$, $i \in [d]$, and $j \in [m]$.

\begin{lemma}{Expected Squared Error}{fast-mmult-variance}
For any sampling probabilities $p_1,\ldots,p_d$, we have
\[
	\E[\normof{\mC-\mA^\intercal\mB}_F^2] \leq \frac1k \sum_{\ell=1}^n \frac1{p_\ell} \normof{\va_\ell}_2^2 \normof{\vb_\ell}_2^2
\]
where $\va_\ell$ and $\vb_\ell$ are the $\ell^{th}$ rows of $\mA$ and $\mB$ respectively.
\end{lemma}
\begin{dropdown}{_Proof_}
\begin{proof}
Let $\mR_t = \frac1{k p_{s_t}} \va_{s_t} \vb_{s_t}^\intercal$, so that we have $\mC = \sum_{t=1}^k \mR_t$:
\begin{align*}
	\mC
	&= (\mS\mA)^\intercal(\mS\mB) \\
	&= \sum_{t=1}^k \frac1{\sqrt{k p_{s_t}}} \va_{s_t} \cdot \frac1{\sqrt{k p_{s_t}}} \vb_{s_t}^\intercal \\
	&= \sum_{t=1}^k \mR_t
\end{align*}
In particular, we see that $\E[\mR_t] = \sum_{\ell=1}^n p_\ell \frac{1}{kp_\ell} \va_\ell\vb_\ell^\intercal = \frac1k \mA^\intercal\mB$, which in turn implies $\E[\mC] = \mA^\intercal\mB$.
We then can expand and simplify by independence, linearity of variance and by the bound $\Var[x] \leq \E[x^2]$:
\begin{align*}
	\E[\normof{\mC-\mA^\intercal\mB}_F^2]
	&= \sum_{i=1}^d\sum_{j=1}^m \E\left[ ([\mC-\mA^\intercal\mB]_{i,j})^2 \right] \\
	&= \sum_{i=1}^d\sum_{j=1}^m \E\left[ \left(\textstyle{\sum_{t=1}^k [\mR_t]_{i,j} - \E[\mR_t]_{i,j}}\right)^2 \right] \\
	&= \sum_{i=1}^d\sum_{j=1}^m \Var\left[ \textstyle{\sum_{t=1}^k [\mR_t]_{i,j}} \right] \\
	&= k \sum_{i=1}^d\sum_{j=1}^m \Var\left[ [\mR_1]_{i,j} \right] \\
	&\leq k \sum_{i=1}^d\sum_{j=1}^m \E\left[ \left([\mR_1]_{i,j}\right)^2 \right] \\
	&\leq k \E\left[ \normof{\mR_1}_F^2 \right]
\end{align*}
Since $\mR_1$ is rank-one, it is simple to compute its Frobenius norm:
\begin{align*}
	\normof{\mR_1}_F^2 
	&= \tr(\mR_1^\intercal\mR_1) \\
	&= \frac1{k^2 p_{s_t}^2} \tr((\va_{s_t} \vb_{s_t}^\intercal)^\intercal(\va_{s_t}\vb_{s_t}^\intercal)) \\
	&= \frac1{k^2 p_{s_t}^2} \tr(\vb_{s_t}\va_{s_t}^\intercal\va_{s_t}\vb_{s_t}^\intercal) \\
	&= \frac{\normof{\va_{s_t}}_2^2 \normof{\vb_{s_t}}_2^2}{k^2 p_{s_t}^2}
\end{align*}
And overall, we conclude that
\begin{align*}
	\E[\normof{\mC-\mA^\intercal\mB}_F^2]
	&\leq k \E\left[ \normof{\mR_1}_F^2 \right] \\
	&= k \cdot \sum_{\ell=1}^n p_\ell \frac{\normof{\va_\ell}_2^2 \normof{\vb_\ell}_2^2}{k^2 p_\ell^2} \\
	&= \frac1k \cdot \sum_{\ell=1}^n \frac{\normof{\va_\ell}_2^2 \normof{\vb_\ell}_2^2}{p_\ell}
\end{align*}
\end{proof}
\end{dropdown}

Having completed this core technical claim, \theoremref{fast-mmult} follows by a short corollary from just plugging in the chosen sampling probabilities.
In fact, we prove something slightly broader:
\begin{corollary}{Oversampling Fast Multiplication}{fast-mmult-oversample}
Let $\tilde\tau_1,\ldots,\tilde\tau_n$ be numbers such that $\tilde\tau_\ell \geq \normof{\va_\ell}_2^2$ for all $\ell\in[n]$, and let $T \defeq \sum_{\ell=1}^n \tilde\tau_\ell$.
Then let $p_\ell \defeq \frac{\tilde\tau_\ell}{T}$ and run \algorithmref{fast-mmult}.
Then, so long as $k \geq \frac{1}{\eps^2 \delta} \cdot \frac{T}{\normof{\mA}_F^2}$, with probability $1-\delta$ we get
\[
	\normof{\mC-\mA^\intercal\mB}_F
	\leq \eps\normof{\mA}_F \normof{\mB}_F
\]
\end{corollary}
\begin{dropdown}{_Proof_}
\begin{proof}
We first bound the expected error from \lemmaref{fast-mmult-variance}, where we get
\begin{align*}
	\normof{\mC-\mA^\intercal\mB}_F^2
	&\leq \frac1k \cdot \sum_{\ell=1}^n \frac{\normof{\va_\ell}_2^2 \normof{\vb_\ell}_2^2}{p_\ell} \\
	&= \frac{T}{k} \cdot \sum_{\ell=1}^n \frac{\normof{\va_\ell}_2^2}{\tilde\tau_\ell} \normof{\vb_\ell}_2^2 \\
	&\leq \frac{T}{k} \cdot \sum_{\ell=1}^n \normof{\vb_\ell}_2^2 \\
	&= \frac{T}{k} \normof{\mB}_F^2
\end{align*}
We apply Markov's inequality, which tells us that
\begin{align*}
	\Pr[\normof{\mC-\mA^\intercal\mB}^2 > \eps^2\normof{\mA}_F^2\normof{\mB}_F^2]
	&\leq \frac{\normof{\mC-\mA^\intercal\mB}_F^2}{\eps^2\normof{\mA}_F^2\normof{\mB}_F^2} \\
	&\leq \frac{\frac 1k T\normof{\mB}_F^2}{\eps^2\normof{\mA}_F^2\normof{\mB}_F^2} \\
	&= \frac{T}{k\eps^2\normof{\mA}_F^2}
\end{align*}
Which is at most $\delta$ when $k > \frac{1}{\delta\eps^2} \cdot \frac{T}{\normof{\mA}_F^2}$.
\end{proof}
\end{dropdown}
Note that when we compute the norms exactly, so that $\tilde\tau_\ell = \normof{\va_\ell}_2^2$ for all $\ell$, we get $T = \sum_\ell \tilde\tau_\ell = \normof{\mA}_F^2$, which recovers \theoremref{fast-mmult} exactly.

# See Also

The analysis here is a blend of \cite{drineas2018lectures} and \cite{nelson15lecture}, as well as some personal notes by Christopher Musco.
Note that randomized matrix multiplication often does not use the exact sampling probabilities discussed here, and the references below discuss a variety of slightly different schemes.

Here's some relevant papers:
- \cite{drineas2006fast} is (afaik) the original paper on the topic. Table 1 of this paper compares a great variety of sampling schemes.
- \cite{drineas2018lectures} is a book with a section that this page partially copies.
- \cite{nelson15lecture} are lecture notes that this page partially copies.
- \cite{avron2019universal} generalizes this to approximate linear operator multiplication in Claim 45, the "Approximate Operator Application".
- _Let me know if anything is missing_

# Bibliography

* \biblabel{avron2019universal}{Avron et al. (2019)} **Avron**, **Kapralov**, **Musco**, **Musco**, **Velingker**, and **Zandieh**. [A Universal Sampling Method for Reconstructing Signals with Simple Fourier Transforms](https://arxiv.org/pdf/1812.08723.pdf). _STOC_ 2019.

* \biblabel{drineas2006fast}{Drineas, Kannan, Mahoney (2006)} **Drineas**, **Kannan**, and **Mahoney**. [Fast Monte Carlo Algorithms for Matrices I: Approximating Matrix Multiplication](https://epubs.siam.org/doi/pdf/10.1137/S0097539704442684). _SIAM Journal on Computing 2006_.

* \biblabel{drineas2018lectures}{Drineas, Mahoney (2018)} **Drineas** and **Mahoney**. [Lectures on Randomized Numerical Linear Algebra](https://arxiv.org/pdf/1712.08880.pdf). _The Mathematics of Data 2018_.

* \biblabel{nelson15lecture}{Nelson (2015)} **Nelson**. [Lecture 15](http://people.seas.harvard.edu/~minilek/cs229r/fall15/lec/lec15.pdf). _Lecture Notes, 2015_.

\theoremscripts
