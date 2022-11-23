@def title = "RandNLA Proof Wiki"
@def tags = ["syntax", "code"]

\enabletheorems
\newcounter{NumAlgorithms}


# Active $L_2$ Regression via Leverage Score Sampling

Consider trying to solve $\min_{\vx} \normof{\mA\vx-\vb}_2$, but under the \highlight{active regression} regime.
That is, we see the full matrix $\mA$ a priori, but do not know the corresponding entries of $\vb$.
We can observe entries of $\vb$, but this is expensive, and we want to read as few entries as possible.

_Prerequisite: [Subspace Embedding via Leverage Score Sampling (Matrix Bernstein)](/leverage-subspace-embedding/)_

_Prerequisite: [Fast Matrix-Matrix Multiplication](/fast-matrix-mult/)_

Without this constraint on the information in $\vb$, solving $\min_{\vx} \normof{\mA\vx-\vb}_2$ from subsampling a small number of rows is relatively straightforward: build a subspace embedding for the rows of the matrix $\bar{\mA} \defeq [\mA \, \vb]$, sampling by the leverage scores of the expanded matrix.

We will analyze the following algorithm:

\begin{algorithm}{Leverage Scores for Active Regression}{leverage-regression}
**input**: Matrix $\mA\in\bbR^{n \times d}$. Query access to $\vb\in\bbR^{n}$. Number $k$ of subsamples.

**output**: Approximate $\tilde{\vx}\in\bbR^{k \times d}$

1. Sample indices $s_1,\ldots,s_k\in[n]$ iid with respect to the leverage scores of $\mA$
1. Let $p_i \defeq \frac{\tau_i}{d}$ denote the probability of sampling row $i$
1. Build the sample-and-rescale matrix $\mS\in\bbR^{k \times n}$:
    
    Row $i$ of $\mS$ has form $\begin{bmatrix}0&0&\cdots&0&\frac{1}{\sqrt{k p_{s_i}}}&0&\cdots&0\end{bmatrix}$, where index $s_i$ is the nonzero entry.

1. Return $\tilde{\vx} = \argmin_{\vx}\normof{\mS\mA-\mS\vb}_2$
\end{algorithm}

Crucially, this algorithm only looks at the rows of $\vb$ selected by the sampling matrix, and this sampling matrix does not know anything about $\vb$.
Hence, the algorithm fits in the active regression model.
We show the following:

\begin{theorem}{Active $L_2$ Regression}{active-l2-regression}
Let $\mA\in\bbR^{n \times d}$ be a tall-and-skinny matrix, and let $\vb\in\bbR^{n}$.
Then let $\tilde{\vx}$ be the result of \algorithmref{leverage-regression} run with $k = \Omega(d \log(\frac d\delta) + \frac{d}{\delta\eps})$.
With probability $1-\delta$, we then have
\[
	\normof{\mA\tilde{\vx}-\vb}_2 \leq (1+\eps) \min_{\vx\in\bbR^{d}}\normof{\mA\vx-\vb}_2
\]
\end{theorem}
\begin{dropdown}{_Proof_}
\begin{proof}
Let $\vx^* \defeq \argmin_{\vx} \normof{\mA\vx-\vb}_2$ be the true solution to the full $L_2$ regression problem.
Recall by setting the derivative of $\normof{\mA\vx-\vb}_2^2$ to be zero, we find that $\mA^\intercal\mA\vx^*-\mA^\intercal\vb=0$.
Equivalently, $\mA^\intercal(\mA\vx^*-\vb)=0$, or in other words, $\mA\vx^*-\vb$ is orthogonal to the range of $\mA$.
In particular, $\mA\vx^*-\vb$ is orthogonal to $\mA\tilde{\vx}-\mA\vx^*$, and so we get
\[
	\normof{\mA\tilde{\vx}-\vb}_2^2 = \normof{\mA\vx^*-\vb}_2^2 + \normof{\mA\tilde{\vx}-\mA\vx^*}_2^2
\]
Let $\mU\in\bbR^{n \times d}$ be an orthogonal matrix than spans the range of $\mA$ (e.g. the $\mU$ matrix from the SVD of $\mA$).
Then, for every $\vx\in\bbR^d$ there exists some $\vy\in\bbR^d$ such that $\mA\vx=\mU\vy$.
Let $\tilde{\vy}$ and $\vy^*$ be the vectors so that $\mA\tilde{\vx}=\mU\tilde{\vy}$ and $\mA\vx^*=\mU\vy^*$.
Then, the above inequality can be written as
\[
	\normof{\mU\tilde{\vy}-\vb}_2^2 = \normof{\mU\vy^*-\vb}_2^2 + \normof{\mU\tilde{\vy}-\mU\vy^*}_2^2
\]
Since $\mU$ is an orthogonal matrix, we have $\normof{\mU\tilde{\vy}-\mU\vy^*}_2^2 = \normof{\tilde{\vy}-\vy^*}_2^2$.
Our goal is to bound the right hand side by $(1+\eps)\normof{\mA\vx^*-\vb}_2^2$, and so it suffices to show that $\normof{\tilde{\vy}-\vy^*}_2^2 \leq \eps \normof{\mU\vy^*-\vb}$.
Recall from [Lemma 1 of the Subspace Embedding Page](/leverage-subspace-embedding/#subspace-equiv) that since $k = \Omega(d \log(\frac d\delta))$, with probability $1-\frac\delta2$, we have $\normof{\mU^\intercal\mS^\intercal\mS\mU-\mI}_2 \leq \frac12$.
Then,
\begin{align*}
	\normof{\tilde{\vy}-\vy^*}_2
	&\leq \normof{(\mU^\intercal\mS^\intercal\mS\mU)(\tilde{\vy}-\vy^*)}_2 + \normof{(\mU^\intercal\mS^\intercal\mS\mU - \mI)(\tilde{\vy}-\vy^*)}_2 \\
	&\leq \normof{(\mU^\intercal\mS^\intercal\mS\mU)(\tilde{\vy}-\vy^*)}_2 + \normof{\mU^\intercal\mS^\intercal\mS\mU - \mI}_2\normof{\tilde{\vy}-\vy^*}_2 \\
	&\leq \normof{(\mU^\intercal\mS^\intercal\mS\mU)(\tilde{\vy}-\vy^*)}_2 + \tsfrac{1}{2} \normof{\tilde{\vy}-\vy^*}_2 \\
	\normof{\tilde{\vy}-\vy^*}_2 &\leq 2\normof{(\mU^\intercal\mS^\intercal\mS\mU)(\tilde{\vy}-\vy^*)}_2
\end{align*}
Next, since $\tilde{\vy} = \argmin_{\vy}\normof{\mS\mU\vy-\mS\vb}_2$, we know that $\mS\mU\tilde{\vy} - \mS\vb$ is orthogonal to everything in the range of $\mS\mU$.
That is $(\mS\mU)^\intercal(\mS\mU\tilde{\vy}-\mS\vb)=0$, so we get that
\begin{align*}
	\normof{(\mU^\intercal\mS^\intercal\mS\mU)(\tilde{\vy}-\vy^*)}_2
	&= \normof{(\mS\mU)^\intercal(\mS\mU\tilde{\vy}-\mS\mU\vy^*)}_2 \\
	&= \normof{(\mS\mU)^\intercal(\mS\mU\tilde{\vy}-\mS\vb+\mS\vb-\mS\mU\vy^*)}_2 \\
	&= \normof{0 + (\mS\mU)^\intercal(\mS\vb-\mS\mU\vy^*)}_2 \\
	&= \normof{(\mS\mU)^\intercal \mS(\mU\vy^*-\vb)}_2
\end{align*}
Where this norm at the end now matches the shape of the [Fast Matrix-Matrix Multiplication](/fast-matrix-mult/) algorithm we studied.
In particular, this is like fast matrix multiplication of the matrices $\mU$ and $\mU\vy^*-\vb$, where the sampling probabilities are proportional to the leverage scores of $\mA$, which by [Lemma 3 of the Basic Properties of Leverage Scores Page](/leverage-score-properties/#basic-props), are equal to the row norms of $\mU$.
So, by [Theorem 1 of the Fast Matrix Multiplication page](fast-matrix-mult/#fast-mmult) with $k = \Omega(\frac{d}{\delta\eps})$ and with probability $1-\frac\delta2$, we have
\[
	\normof{(\mS\mU)^\intercal \mS(\mU\vy^*-\vb)}_2 \leq \frac{\sqrt{\eps}}{2\sqrt{d}} \normof{\mU}_F \normof{\mU\vy^*-\vb}_2
\]
$\mU$ is an orthogonal matrix, so each column of $\mU$ has norm 1, so $\normof{\mU}_F^2 = d$.
Therefore, we get
\[
	\normof{(\mS\mU)^\intercal \mS(\mU\vy^*-\vb)}_2 \leq \frac{\sqrt{\eps}}{2} \normof{\mU\vy^*-\vb}_2
\]
Backing up, we have overall proven that
\begin{align*}
	\normof{\mA\tilde{\vx}-\vb}_2^2
	&= \normof{\mA\vx^*-\vb}_2^2 + \normof{\mA\tilde{\vx}-\mA\vx^*}_2^2 \\
	&= \normof{\mU\vy^*-\vb}_2^2 + \normof{\mU\tilde{\vy}-\mU\vy^*}_2^2 \\
	&= \normof{\mU\vy^*-\vb}_2^2 + \normof{\tilde{\vy}-\vy^*}_2^2 \\
	&\leq \normof{\mU\vy^*-\vb}_2^2 + 4\normof{(\mU^\intercal\mS^\intercal\mS\mU)(\tilde{\vy}-\vy^*)}_2^2 \\
	&\leq \normof{\mU\vy^*-\vb}_2^2 + \eps\normof{\mU\vy^*-\vb}_2^2 \\
	&= (1+\eps)\normof{\mU\vy^*-\vb}_2^2 \\
	&= (1+\eps)\normof{\mA\vx^*-\vb}_2^2
\end{align*}
Which completes the proof.
\end{proof}
\end{dropdown}

---

Notably, this proof uses the sampling probabilities in only two places: showing the subspace embedding guarantee $\normof{\mU^\intercal\mS^\intercal\mS\mU-\mI}_2\leq\frac12$ and showing the fast matrix multiplication $\normof{(\mS\mU)^\intercal \mS(\mU\vy^*-\vb)}_2 \leq \frac{\sqrt{\eps}}{2\sqrt{d}} \normof{\mU}_F \normof{\mU\vy^*-\vb}_2$.
Both of these proofs work with inexact leverage scores sampling probabilities, namely via oversampling, as shown on their respective pages ([Theorem 3 here](/leverage-subspace-embedding/#oversampling) and [Corollary 1 here](/fast-matrix-mult/#fast-mmult-oversample)).
We can therefore state a more general algorithm and guarantee for leverage score sampling with the exact same proof analysis:

\begin{algorithm}{Leverage Scores for Active Regression}{leverage-regression-oversampling}
**input**: Matrix $\mA\in\bbR^{n \times d}$. Query access to $\vb\in\bbR^{n}$. Number $k$ of subsamples. Sampling probabilities $p_1,\ldots,p_n$.

**output**: Approximate $\tilde{\vx}\in\bbR^{k \times d}$

1. Sample indices $s_1,\ldots,s_k\in[n]$ iid with respect to $p_1,\ldots,p_n$
1. Let $p_i \defeq \frac{\tau_i}{d}$ denote the probability of sampling row $i$
1. Build the sample-and-rescale matrix $\mS\in\bbR^{k \times n}$:
    
    Row $i$ of $\mS$ has form $\begin{bmatrix}0&0&\cdots&0&\frac{1}{\sqrt{k p_{s_i}}}&0&\cdots&0\end{bmatrix}$, where index $s_i$ is the nonzero entry.

1. Return $\tilde{\vx} = \argmin_{\vx}\normof{\mS\mA-\mS\vb}_2$
\end{algorithm}

\begin{theorem}{Active $L_2$ Regression (Oversampling)}{active-l2-regression-oversampling}
Fix $\eps>0,\delta>0$.
Let $\tau_1,\ldots,\tau_n$ be the leverage scores of $\mA\in\bbR^{n \times d}$, and let $\tilde\tau_1,\ldots,\tilde\tau_n$ be _overestimates_ of the leverage scores, so that $\tilde\tau_i \geq \tau_i$.
Let $\tilde d \defeq \sum_{i=1}^n \tilde\tau_i$.
Then let $\tilde{\vx}$ be the output of \algorithmref{active-l2-regression-oversampling} with $k = \Omega(\tilde d \log(\frac d \delta) + \frac{\tilde d}{\eps \delta})$ and probabilities $p_i = \frac{\tilde\tau_i}{\tilde d}$
With probability $1-\delta$, we then have
\[
	\normof{\mA\tilde{\vx}-\vb}_2^2 \leq (1+\eps) \min_{\vx\in\bbR^{d}}\normof{\mA\vx-\vb}_2^2
\]
\end{theorem}
\begin{dropdown}{_Proof_}
\begin{proof}
To get the subspace embedding, we directly appeal to [Theorem 3 here](/leverage-subspace-embedding/#oversampling).
We trace through the fast matrix multiplication more carefully.
We appeal to [Corollary 1 here](/fast-matrix-mult/#fast-mmult-oversample) with error tolerance $\eps_0 = \sqrt{\frac\eps d}$, and with $T = \tilde d$ since $\tilde\tau_i$ overestimate the row norms of $\mU$.
So, we need to sample
\[
	k
	\geq \frac{1}{\eps_0^2 \delta} \cdot \frac{T}{\normof{\mU}_F^2}
	= \frac{d}{\eps \delta} \cdot \frac{\tilde d}{d}
	= \frac{\tilde d}{\eps \delta}
\]
rows to maintain correctness.
\end{proof}
\end{dropdown}

# See Also

The proofs above are a near-exact copy of an unpublished note by Christopher Musco, and it closely tracks an analysis of \cite{woodruff2014sketching} which studies regression under oblivious sketching instead of leverage score sampling.

There are resources with similar analyses in the literature:
- \citep{sarlos2006improved} Is the original paper and analysis.
- \citep{chen2019active} Explicitly studies active regression under Leverage Functions, which relates to Linear Operators, more generally than just matrices.
- \citep{musco2021active} Studies active regression beyond the $\ell_2$ norm, like in the $\ell_p$ norm.
- \citep{meyer2020statistical} The appendix studies active $\ell_2$ regression for both matrices and operators, but where we choose one of many design matrices $\mA_1,\ldots,\mA_q$, and can only achieve constant factor error (i.e. $\eps > \Omega(1)$).
- \citep{kapralov2023toeplitz} Contains an analysis of active $\ell_2$ regression, but where we choose one of many design matrices $\mA_1,\ldots,\mA_q$, and can achieve $(1+\eps)$ factor error, but with a $\tilde O(\tilde d^2)$ dependence.
- _Let me know if anything is missing_

# Bibliography

* \biblabel{woodruff2014sketching}{Woodruff (2014)} **Woodruff**. [Sketching as a Tool for Numerical Linear Algebra.](http://www.cs.cmu.edu/afs/cs/user/dwoodruf/www/wNow3.pdf). _Foundations and Trends® in Theoretical Computer Science_ 2014.

* \biblabel{sarlos2006improved}{Sarlos (2017)} **Sarlos**. [Improved Approximation Algorithms for Large Matrices via Random Projections](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4031351). _FOCS_ 2006.

* \biblabel{chen2019active}{Chen Price (2019)} **Chen** and **Price**. [Active Regression via Linear-Sample Sparsification](https://arxiv.org/pdf/1511.07263.pdf). _COLT_ 2019.

* \biblabel{musco2021active}{Musco et al. (2019)} **Musco**, **Musco**, **Woodruff**, and **Yasuda**. [Active Sampling for Linear Regression Beyond the $\ell_2 $ Norm](https://arxiv.org/pdf/2111.04888.pdf). _FOCS_ 2022.

* \biblabel{meyer2020statistical}{Musco et al. (2019)} **Meyer** and **Musco**. [The Statistical Cost of Robust Kernel Hyperparameter Turning](https://arxiv.org/pdf/2006.08035.pdf). _NeurIPS_ 2020.

* \biblabel{kapralov2023toeplitz}{Kapralov et al. (2023)} **Kapralov**, **Lawrence**, **Makarov**, **Musco**, and **Sheth**. [Toeplitz Low-Rank Approximation with Sublinear Query Complexity](https://arxiv.org/pdf/2211.11328.pdf). _SODA_ 2023.



\theoremscripts
