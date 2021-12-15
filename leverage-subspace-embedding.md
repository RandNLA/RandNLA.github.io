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

\begin{lemma}{Equivalent form of Subspace Embedding}{subspace-equiv}
Let $\mU\in\bbR^{n \times d}$ be any matrix with orthonormal columns that space the columnspace of $\mA$.
Then the condition in \theoremref{leverage-subspace-embedding} is equivalent to
\[
	\onormof{\mU^\intercal\mS^\intercal\mS\mU - \mI}{2}\leq\eps
\]
\end{lemma}
\begin{dropdown}{_Proof_}
\begin{proof}
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
\end{proof}
\end{dropdown}
Note that this guarantee seems reasonable.
Below, we verify that $\E[\mS^\intercal\mS]=\mI$, so we have $\E[\mU^\intercal\mS^\intercal\mU\mS] = \mU^\intercal\mU = \mI$, and thus $\normof{\E[\mU^\intercal\mS^\intercal\mS\mU] - \mI}_2 = 0$.
In other words, this guarantee looks like a concentration of a random variable around its mean.


\begin{lemma}{Unbiased Expectation}{expectation}
The Sample and Rescale matrix $\mS$ has $\E[\mS^\intercal\mS]=\mI$.
\end{lemma}
\begin{dropdown}{_Proof_}
\begin{proof}
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
\end{proof}
\end{dropdown}

---

With this reformulation, we now turn to ensuring the Subspace Embedding guarantee $\normof{\mU^\intercal\mS^\intercal\mS\mU - \mI}_2 \leq \eps$
In particular, we import the following theorem:

\begin{theorem}{Matrix Bernstein (Simplified)}{matrix-bernstein}
Let $\mR_1,\ldots,\mR_k\in\bbR^{d \times d}$ be iid symmetric random matrices such that $\E[\mR]=0$, $\normof{\mR}_2\leq\gamma$, and $\normof{\E[\mR^2]}_2\leq\sigma^2$.
Then,
\[
	\Pr[\normof{\frac1k\sum_{i=1}^k \mR_i}_2 \geq \eps] \leq 2de^{\frac{-k\eps^2}{2\sigma^2+\frac{2}{3}\gamma\eps}}
\]
_(This is a heavy simplification of Theorem 6.1.1 from \cite{tropp2015introduction})_
\end{theorem}

We now prove the main theorem:
\begin{dropdown}{_Proof of \theoremref{leverage-subspace-embedding}_}
\begin{proof}
We now decompose $\mI-\mU^\intercal\mS^\intercal\mS\mU = \frac1k\sum_{i=1}^k \mR_i$ for some matrices $\mR_i$.
First, note that the $i^{th}$ row of $\mS\mU$ is $\frac{1}{\sqrt{k p_{s_i}}} \vu_{s_i}$, where $\vu_j$ denotes the $j^{th}$ row of $\mU$.
In particular, by the outer-product view of matrix-matrix products, note that
\[
	\mU^\intercal\mS^\intercal\mS\mU = (\mS\mU)^\intercal(\mS\mU) = \sum_{i=1}^k \tfrac{1}{kp_{s_i}} \vu_{s_i}\vu_{s_i}^\intercal
\]
Which then motivates us to let $\mR_i \defeq \mI - \frac{1}{p_{s_i}} \vu_{s_i}\vu_{s_i}^\intercal$ (notice $k$ is dropped from the denominator) be our iid matrices.
- By the above equation, we know that $\frac1k\sum_{i=1}^k\mR_i = \mI-\mU^\intercal\mS^\intercal\mS\mU$.
- $\E[\mR] = \mI - \sum_{j=1}^n p_j \frac1{p_j}\vu_j\vu_j^\intercal = \mI - \sum_{j=1}^n \vu_j\vu_j^\intercal = \mat{0}$
- For all $j\in[n]$,  $\normof{\mR}_2 \leq \normof{\mI}_2 + \frac1{p_j}\normof{\vu_j\vu_j^\intercal}_2 = 1 + \frac{d}{\normof{\vu_j}_2^2} \normof{\vu_j}_2^2 = 1 + d$ (using [Lemma 3 Part 4 from here](/leverage-score-properties/#basic-props))
- We verify the variance term $\normof{\E[\mR^2]}_2 = (d-1)\mI$ below:
\begin{align}
	\E[\mR^2]
	&= \E[\mI - \tfrac2{p_j} \vu_j\vu_j^\intercal + \tfrac1{p_j^2} \vu_j\vu_j^\intercal\vu_j\vu_j^\intercal] \\
	&= \mI - 2\mI + \E[\tfrac1{p_j^2} \vu_j\vu_j^\intercal\vu_j\vu_j^\intercal] \\
	&= \mI - 2\mI + \sum_{j=1}^n[p_j \tfrac1{p_j^2} \vu_j\vu_j^\intercal\vu_j\vu_j^\intercal] \\
	&= \mI - 2\mI + \sum_{j=1}^n[\tfrac1{p_j} \normof{\vu_j}_2^2 \vu_j\vu_j^\intercal] \\
	&= \mI - 2\mI + \sum_{j=1}^n[\tfrac{d}{\normof{\vu_j}_2^2} \normof{\vu_j}_2^2 \vu_j\vu_j^\intercal] \\
	&= \mI - 2\mI + d\sum_{j=1}^n[\vu_j\vu_j^\intercal] \\
	&= \mI - 2\mI + d\mI \\
	&= (d-1)\mI
\end{align}

And so, by the Matrix Bernstein bound,
\[
	\Pr[\normof{\mI-\mU^\intercal\mS^\intercal\mS\mU}_2 \geq \eps] \leq 2d e^{-\frac{k\eps^2}{2(d-1) + \frac23 (d+1)\eps}}
\]
Plugging in $k = O(\frac{d \log(\frac d\delta)}{\eps^2})$ makes the right hand side term at most $\delta$, completing the proof.
\end{proof}
\end{dropdown}

# Oversampling: Using approximate leverage scores

A very important and common extension of row sampling is discussed here.
Note that \theoremref{leverage-subspace-embedding} only works, as stated, if we compute the leverage scores of $\mA$ exactly.
This is not necessary, and approximate leverage scores suffice, as characterized below:

\begin{algorithm}{Leverage Score Subspace Embedding}
**input**: Matrix $\mA\in\bbR^{n \times d}$. Number $k$ of subsamples. Sampling probabilities $p_1,\ldots,p_n$.

**output**: Sketched matrix $\tilde{\mA}\in\bbR^{k \times d}$

1. Sample indices $s_1,\ldots,s_k\in[n]$ iid with respect to $p_i,\ldots,p_n$
1. Build the sample-and-rescale matrix $\mS\in\bbR^{k \times n}$:
    
    Row $i$ of $\mS$ has form $\begin{bmatrix}0&0&\cdots&0&\frac{1}{\sqrt{k p_{s_i}}}&0&\cdots&0\end{bmatrix}$, where index $s_i$ is the nonzero entry.

1. Return $\tilde{\mA} = \mS\mA$
\end{algorithm}

\nbsp

\begin{theorem}{Subspace Embedding (Oversampling)}{oversampling}
Fix $\eps > 0$, $\delta > 0$.
Let $\tau_1,\ldots,\tau_n$ be the leverage scores for $\mA\in\bbR^{n \times d}$, and let $\tilde\tau_1,\ldots,\tilde\tau_n$ be _overestimates_ of the leverage scores, so that $\tilde\tau_i \geq \tau_i$.
Let $\tilde d \defeq \sum_{i=1}^n \tilde\tau_i$.
Then run Algorithm 2 with $k = O(\frac{\tilde d\log(\frac d\delta)}{\eps^2})$ and sampling probabilities $p_i = \frac{\tilde\tau_i}{\tilde d}$.
Then with probability $1-\delta$, we have _for all_ $\vx\in\bbR^{d}$
\[
	(1-\eps)\normof{\mA\vx}_2^2 \leq \normof{\mS\mA\vx}_2^2 \leq (1+\eps)\normof{\mA\vx}_2^2
\]
\end{theorem}
\begin{dropdown}{_Proof_}
\begin{proof}
This proof almost exactly matches that of \theoremref{leverage-subspace-embedding}, except we bound a couple terms from the Matrix Bernstein a little differently:
- $\normof{\mR}_2 \leq \normof{\mI}_2 + \frac1{p_j}\normof{\vu_j\vu_j^\intercal}_2 = 1 + \frac{\tilde d}{\tilde\tau_i} \normof{\vu_j}_2^2 = 1 + \frac{\tilde d}{\tilde \tau_i} \tau_i \leq 1 + \tilde d$
- And the variance analysis:
\begin{align}
	\E[\mR^2]
	&= \E[\mI - \tfrac2{p_j} \vu_j\vu_j^\intercal + \tfrac1{p_j^2} \vu_j\vu_j^\intercal\vu_j\vu_j^\intercal] \\
	&= \mI - 2\mI + \E[\tfrac1{p_j^2} \vu_j\vu_j^\intercal\vu_j\vu_j^\intercal] \\
	&= \mI - 2\mI + \sum_{j=1}^n[p_j \tfrac1{p_j^2} \vu_j\vu_j^\intercal\vu_j\vu_j^\intercal] \\
	&= \mI - 2\mI + \sum_{j=1}^n[\tfrac1{p_j} \normof{\vu_j}_2^2 \vu_j\vu_j^\intercal] \\
	&= \mI - 2\mI + \sum_{j=1}^n[\tfrac{\tilde d}{\tilde\tau_i} \tau_i \vu_j\vu_j^\intercal] \\
	&\leq \mI - 2\mI + \tilde d\sum_{j=1}^n[\vu_j\vu_j^\intercal] \\
	&= \mI - 2\mI + \tilde d\mI \\
	&= (\tilde d-1)\mI
\end{align}
From which the result follows by Matrix Bernstein.
\end{proof}
\end{dropdown}

Recall that the true leverage scores have $\sum_{i=1}^n\tau_i = d$, so $\tilde d$ in \theoremref{oversampling} measures how inaccurate our approximate leverage scores are.
A nice example to pair with this theorem is how it explain when _uniform_ sampling works:

\begin{corollary}{Subspace Embedding via Uniform Sampling}{uniform-sampling}
Let $\tau^* \defeq \max_{i\in[n]} \tau_i(\mA)$.
Then run Algorithm 2 with $k = O(\frac{n\tau^* \log(\frac{d}{\delta})}{\eps^2})$ by with uniform probabilities $p_i = \frac1n$.
Then, with probability $1-\delta$, we have _for all_ $\vx\in\bbR^{d}$
\[
	(1-\eps)\normof{\mA\vx}_2^2 \leq \normof{\mS\mA\vx}_2^2 \leq (1+\eps)\normof{\mA\vx}_2^2
\]
\end{corollary}

Matrices with large $\tau^*\approx1$ have a unique row that is not well represented by any other row in $\mA$, so $\tilde O(n)$ row samples are necessary to sample that one unique row.
Matrices with low $\tau^*\approx\frac1n$ have rows that all looks almost equivalent, so just sampling $\tilde O(1)$ rows suffices to effectively recover all of $\mA$.
This explain when uniform sampling can be better or worse than leverage score sampling, depending on the exact matrix $\mA$.
This is also why _(in)coherence_ assumptions often appear in the Matrix Completion literature.

# See Also

There are many extensions and implementations of subspace embedding in the literature:
- _Let me know if anything is missing_

# Bibliography

* \biblabel{tropp2015introduction}{Tropp (2015)} **Tropp**. [An Introduction to Matrix Concentration Inequalities](http://users.cms.caltech.edu/~jtropp/books/Tro14-Introduction-Matrix-FnTML-rev.pdf). _Foundations and Trends in Machine Learning_ 2015.


\theoremscripts
