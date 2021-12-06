@def title = "RandNLA Proof Wiki"
@def tags = ["syntax", "code"]

\newcounter{NumAlgorithms}

\enabletheorems


# Subspace Embedding via Leverage Score Sampling
On this page, we show how sampling a matrix via its leverage scores creates a subspace embedding guarantee.
Subspace embeddings are core to much of RandNLA, so many papers either assume they are given a subspace embedding, or create a proof much like the one below (just fine-tuned to research their setting).

_Prerequisite: [Basic Properties of Leverage Scores](/leverage-score-properties/)_

We consider the following algorithm:
\begin{algorithm}{Leverage Score Subspace Embedding}
**input**: Matrix $\mA\in\bbR^{n \times d}$. Number $k$ of subsamples.

**output**: Sketched matrix $\tilde{\mA}\in\bbR^{k \times d}$

1. Sample indices $s_1,\ldots,s_k\in[n]$ iid with respect to the leverage scores of $\mA$
1. Let $p_i \defeq \frac{\tau_i}{d}$ denote the probability of sampling row $i$
1. Build the sample-and-rescale matrix $\mS\in\bbR^{k \times n}$:
    
    Row $i$ of $\mS$ has form $\begin{bmatrix}0&0&\cdots&0&\frac{1}{\sqrt{k p_{s_i}}}&0&\cdots&0\end{bmatrix}$, where index $s_i$ is the nonzero entry.

1. Return $\tilde{\mA} = \mS\mA$
\end{algorithm}

This algorithm is equivalent to building $\tilde{\mA}$ by sampling $k$ rows of $\mA$ with respect to leverage scores, rescaling those rows by $\frac{1}{\sqrt{k p_{s_i}}}$, and stacking those rows up.
The matrix $\mS$ is just a formalization that makes it easier to analyze this algorithm.

\begin{theorem}{Subspace Embedding}{leverage-subspace-embedding}
Fix $\eps > 0$, $\delta > 0$ and run the above algorithm with $k = O(\frac{d\log(\frac d\delta)}{\eps^2})$.
Then with probability $1-\delta$, we have _for all_ $\vx\in\bbR^{d}$
\[
	(1-\eps)\normof{\mA\vx}_2^2 \leq \normof{\mS\mA\vx}_2^2 \leq (1+\eps)\normof{\mA\vx}_2^2
\]
\end{theorem}

---

To prove this result, we will first massage the _Subspace Embedding Guarantee_ into a more approachable form.
First, let $\mU\in\bbR^{n \times d}$ be any matrix with orthonormal columns with the same columnspace as $\mA$ (e.g. the $\mU$ matrix in the SVD of $\mA$).
Then, for any $\vx\in\bbR^d$ we have a one-to-one map to a vector $\vy\in\bbR^d$ such that $\mA\vx=\mU\vy$.
So, we need to show that
\begin{align}
	(1-\eps)\normof{\mU\vy}_2^2 \leq \normof{\mS\mU\vy}_2^2 \leq (1+\eps)\normof{\mU\vy}_2^2 && \forall\vy\in\bbR^d
\end{align}
Next, we expand the $\ell_2$ norms:
\begin{align}
	(1-\eps)\vy^\intercal\mU^\intercal\mU\vy \leq \vy^\intercal\mU^\intercal\mS^\intercal\mS\mU\vy \leq (1+\eps)\vy^\intercal\mU^\intercal\mU\vy && \forall\vy\in\bbR^d
\end{align}
Note that $\mU^\intercal\mU=\mI$ by definition, so we get
\begin{align}
	(1-\eps)\vy^\intercal\vy \leq \vy^\intercal\mU^\intercal\mS^\intercal\mS\mU\vy \leq (1+\eps)\vy^\intercal\vy && \forall\vy\in\bbR^d
\end{align}
For $\vy=0$, the claim trivially holds, so assuming $\vy\neq0$ we divide by $\vy^\intercal\vy$ on both sides:
\begin{align}
	(1-\eps) \leq \frac{\vy^\intercal\mU^\intercal\mS^\intercal\mS\mU\vy}{\vy^\intercal\vy} \leq (1+\eps) && \forall\vy\in\bbR^d
\end{align}
By the [Courant–Fischer-Weyl Min-Max Principle](https://en.wikipedia.org/wiki/Min-max_theorem), this is equivalent to guaranteeing that
\begin{align}
	1-\eps \leq \lambda_i(\mU^\intercal\mS^\intercal\mS\mU) \leq 1+\eps && \forall i\in[d]
\end{align}
And rearranging the claim, this becomes
\begin{align}
	\abs{\lambda_i(\mU^\intercal\mS^\intercal\mS\mU) - 1} \leq \eps && \forall i\in[d]
\end{align}
Or, equivalently,
\[
	\normof{\mU^\intercal\mS^\intercal\mS\mU - \mI}_2 \leq \eps
\]
Note that this guarantee seems reasonable.
We can verify that $\E[\mS^\intercal\mS]=\mI$, so we have $\E[\mU^\intercal\mS^\intercal\mU\mS] = \mU^\intercal\mU = \mI$, and thus $\normof{\E[\mU^\intercal\mS^\intercal\mS\mU] - \mI}_2 = 0$.
In other words, this guarantee looks like a concentration of a random variable around its mean.

---

We now verify this claim that $\E[\mS^\intercal\mS]=\mI$.
First, we let $\vr_i$ denote row $i$ of $\mS$:
\[
	\vr_i = \begin{bmatrix}0&0&\cdots&0&\frac{1}{\sqrt{k p_{s_i}}}&0&\cdots&0\end{bmatrix}
\]
We can equivalently define $\vr$ entrywise, with an indicator variable: $[\vr_i]_j = \frac1{\sqrt{k p_j}} \indicate{s_i = j}$.
By the outer-product form of matrix-matrix multiplication, note that
\[
	\E[\mS^\intercal\mS] = \E\left[\sum_{i=1}^k \vr_i\vr_i^\intercal \right] = k \E[\vr_1\vr_1^\intercal]
\]
So, it suffices to analyze the expected value of $\vr_1\vr_1^\intercal$.
In particular, note that $\vr_1\vr_1^\intercal$ is a diagonal matrix with entries $\frac1{k p_j} \indicate{s_i = j}$.
This lets us compute the expected value of an arbitrary diagonal entry of $\vr_1\vr_1^\intercal$ then:
\[
	\E[[\vr_1\vr_1^\intercal]_{j,j}]
	= \E[\frac1{k p_j} \indicate{s_i = j}]
	= \frac1{k p_j} \Pr[s_i = j]
	= \frac1{k p_j} p_j
	= \frac1{k}
\]
And so, we have $\E[\vr_1\vr_1^\intercal] = \frac1k\mI$, and thus $\E[\mS^\intercal\mS] = k\E[\vr_1\vr_1^\intercal]=\mI$.

---

With this reformulation, we now turn to ensuring the Subspace Embedding guarantee $\normof{\mU^\intercal\mS^\intercal\mS\mU - \mI}_2 \leq \eps$
In particular, we import the following theorem:

\begin{theorem}{Matrix Chernoff (Simplified)}{matrix-chernoff}
Let $\mR_1,\ldots,\mR_k\in\bbR^{d \times d}$ be iid symmetric random matrices such that $\E[\mR]=0$, $\normof{\mR}_2\leq\gamma$, and $\normof{\E[\mR^2]}_2\leq\sigma^2$.
Then,
\[
	\Pr[\normof{\frac1k\sum_{i=1}^k \mR_i}_2 \geq \eps] \leq 2de^{\frac{-k\eps^2}{2\sigma^2+\frac{2}{3}\gamma\eps}}
\]
\end{theorem}
This is theorem XXX from XXX.

We now decompose $\mI-\mU^\intercal\mS^\intercal\mS\mU = \frac1k\sum_{i=1}^k \mR_i$ for some matrices $\mR_i$.
First, note that the $i^{th}$ row of $\mS\mU$ is $\frac{1}{\sqrt{k p_{s_i}}} \vu_{s_i}$, where $\vu_j$ denotes the $j^{th}$ row of $\mU$.
In particular, by the outer-product view of matrix-matrix products, note that
\[
	\mU^\intercal\mS^\intercal\mS\mU = (\mS\mU)^\intercal(\mS\mU) = \sum_{i=1}^k \tfrac{1}{kp_{s_i}} \vu_{s_i}\vu_{s_i}^\intercal
\]
Which then motivates us to let $\mR_i \defeq \mI - \frac{1}{p_{s_i}} \vu_{s_i}\vu_{s_i}^\intercal$ (notice $k$ is dropped from the denominator).
- By the above equation, we know that $\frac1k\sum_{i=1}^k\mR_i = \mI-\mU^\intercal\mS^\intercal\mS\mU$.
- $\E[\mR_i] = \mI - \sum_{j=1}^n p_j \frac1{p_j}\vu_j\vu_j^\intercal = \mI - \sum_{j=1}^n \vu_j\vu_j^\intercal = \mat{0}$
- $\normof{\mR_i}_2 \leq \normof{\mI}_2 + \frac1{p_j}\normof{\vu_{s_i}\vu_{s_i}^\intercal}_2 = 1 + \frac{d}{\normof{\vu_{s_j}}_2^2} \normof{\vu_{s_j}}_2^2 = 1 + d$ (using [Lemma X from Y](/index/))
- We verify the variance term $\normof{\E[\mR_i^2]}_2 = (d-1)\mI$ below:
\begin{align}
	\E[\mR_i^2]
	&= \E[\mI - \tfrac2{p_{s_i}} \vu_{s_i}\vu_{s_i}^\intercal + \tfrac1{p_{s_i}^2} \vu_{s_i}\vu_{s_i}^\intercal\vu_{s_i}\vu_{s_i}^\intercal] \\
	&= \mI - 2\mI + \E[\tfrac1{p_{s_i}^2} \vu_{s_i}\vu_{s_i}^\intercal\vu_{s_i}\vu_{s_i}^\intercal] \\
	&= \mI - 2\mI + \sum_{j=1}^n[p_j \tfrac1{p_j^2} \vu_j\vu_j^\intercal\vu_j\vu_j^\intercal] \\
	&= \mI - 2\mI + \sum_{j=1}^n[\tfrac1{p_j} \normof{\vu_j}_2^2 \vu_j\vu_j^\intercal] \\
	&= \mI - 2\mI + \sum_{j=1}^n[\tfrac{d}{\normof{\vu_j}_2^2} \normof{\vu_j}_2^2 \vu_j\vu_j^\intercal] \\
	&= \mI - 2\mI + d\sum_{j=1}^n[\vu_j\vu_j^\intercal] \\
	&= \mI - 2\mI + d\mI \\
	&= (d-1)\mI
\end{align}

And so, by the matrix Chernoff bound,
\[
	\Pr[\normof{\mI-\mU^\intercal\mS^\intercal\mS\mU}_2 \geq \eps] \leq 2d e^{-\frac{k\eps^2}{2(d-1) + \frac23 (d+1)\eps}}
\]
Plugging in $k = O(\frac{d \log(\frac d\delta)}{\eps^2})$ makes the right hand side term at most $\delta$, completing the proof.

# Oversampling: A core technique for approximate leverage scores

# See Also

The proofs above are most closely related to their highly generalized analysis in Theorem 5 from \cite{avron2019universal}.

There are many extensions of leverage scores in the literature:
- \citep{avron2011input} Ridge Leverage Scores, which relate to Ridge Regression.
- \citep{alaoui2014fast} Kernel Ridge Leverage Scores, which relate to Kernel Ridge Regression.
- \citep{avron2011randomized} Kernel Ridge Leverage Function, which relates to Random Fourier Features for Kernel Ridge Regression.
- \citep{chen2019active} Leverage Function, which relates to Linear Operators instead of matrices. Sometimes called a _Sensitvitiy Score_.
- \citep{avron2019universal} Ridge Leverage Function, which relates to Linear Operators instead of matrices for Ridge Regression .
- \citep{cohen2015lp} Lewis Weights, which relate to $L_p$ norms instead of just $L_2$ norms.
- _Let me know if anything is missing_

# Bibliography

* \biblabel{avron2019universal}{Avron et al. (2019)} **Avron**, **Kapralov**, **Musco**, **Musco**, **Velingker**, and **Zandieh**. [A Universal Sampling Method for Reconstructing Signals with Simple Fourier Transforms](https://arxiv.org/pdf/1812.08723.pdf). _STOC_ 2019.

* \biblabel{alaoui2014fast}{Alaoui Mahoney (2015)} **Alaoui** and **Mahoney**. [Fast Randomized Kernel Methods with Statistical Guarantees.](https://www.stat.berkeley.edu/~mmahoney/pubs/elalaoui-nips15.pdf). _NIPS_ 2015.

* \biblabel{avron2011randomized}{Avron et al. (2017)} **Avron**, **Kapralov**, **Musco**, **Musco**, **Velingker**, and **Zandieh**. [Random Fourier Features for Kernel Ridge Regression: Approximation Bounds and Statistical Guarantees](https://arxiv.org/pdf/1804.09893.pdf). _ICML_ 2017.

* \biblabel{avron2011input}{Cohen Musco² (2017)} **Cohen**, **Kapralov**, **Musco**, and **Musco**. [Input Sparsity Time Low-Rank Approximation via Ridge Leverage Score Sampling](https://arxiv.org/pdf/1511.07263.pdf). _SODA_ 2017.

* \biblabel{chen2019active}{Chen Price (2019)} **Chen** and **Price**. [Active Regression via Linear-Sample Sparsification](https://arxiv.org/pdf/1511.07263.pdf). _COLT_ 2019.

* \biblabel{cohen2015lp}{Cohen Peng (2019)} **Cohen** and **Peng**. [$\ell_p$ Row Sampling by Lewis Weights](https://arxiv.org/pdf/1412.0588.pdf). _STOC_ 2015.


\theoremscripts
