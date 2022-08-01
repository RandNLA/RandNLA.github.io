@def title = "RandNLA Proof Wiki"
@def tags = ["syntax", "code"]

\enabletheorems


# Basic Properties of Leverage Scores
On this page, we will show that three different definitions of Leverage Scores are equivalent.
As a small supplement, we also use two of these forms to bound the maximum value of a leverage score, and also to compute the sum of all leverage scores.

\begin{definition}{Max Characterization}{max-char}
Let $\mA\in\bbR^{n \times d}$.
Then the Maximum Characterization of the Leverage Score for row $i$ of $\mA$ is
\[
	\tau_i(\mA) \defeq \max_{\vx\in\bbR^d} \frac{[\mA\vx]_i^2}{\normof{\mA\vx}_2^2}
\]

\end{definition}

---

\begin{lemma}{Inner Product Characterization}{inner-product-char}
Let $\va_i$ be the $i^{th}$ row of $\mA$.
Then,
\[
	\tau_i(\mA) = \va_i^\intercal(\mA^\intercal\mA)^{-1}\va_i
\]
\end{lemma}
\begin{dropdown}{_Proof_}
\begin{proof}
This is proven in two steps.
First, we show that the inner product characterization upper bounds the max characterization.
Then, we show a matching lower bound.
For simplicity, we prove this for full-rank $\mA$.

Before getting started, we note a useful equation:
\begin{align}
	\normof{\mA(\mA^\intercal\mA)^{-1}\va_i}_2^2
	&= (\va_i^\intercal(\mA^\intercal\mA)^{-1}\mA^\intercal\mA(\mA^\intercal\mA)^{-1}\va_i) \\
	&= (\va_i^\intercal(\mA^\intercal\mA)^{-1}\va_i)
\end{align}

_Upper Bound:_
To create the upper bound, we relate
\begin{align}
	[\mA\vx]_i^2
	&= (\va_i^\intercal\vx)^2 \\
	&= (\va_i^\intercal(\mA^\intercal\mA)^{-1}(\mA^\intercal\mA)\vx)^2 \\
	&= ((\mA(\mA^\intercal\mA)^{-1}\va_i)^\intercal ~ (\mA\vx))^2 \\
	&\leq \normof{\mA(\mA^\intercal\mA)^{-1}\va_i}_2^2 \cdot \normof{\mA\vx}_2^2 \\
	&= (\va_i^\intercal(\mA^\intercal\mA)^{-1}\va_i) \cdot \normof{\mA\vx}_2^2
\end{align}
Where the inequality used is the Cauchy-Schwarz Inequality: $(\vv^\intercal\vy)^2\leq \normof{\vv}_2^2 \cdot \normof{\vy}_2^2$.
We can then give an upper bound to the max characterization:
\[
	\tau_i(\mA)
	= \max_{\vx\in\bbR^d} \frac{[\mA\vx]_i^2}{\normof{\mA}_2^2}
	\leq \max_{\vx\in\bbR^d} \frac{(\va_i(\mA^\intercal\mA)^{-1}\va_i) \cdot \normof{\mA\vx}_2^2}{\normof{\mA\vx}_2^2}
	= \va_i(\mA^\intercal\mA)^{-1}\va_i
\]
_Lower Bound:_
For the lower bound, we just plug $\vx=(\mA^\intercal\mA)^{-1}\va_i$ into the max characterization:
\begin{align}
	\tau_i(\mA)
	\geq \frac{[\mA(\mA^\intercal\mA)^{-1}\va_i]_i^2}{\normof{\mA(\mA^\intercal\mA)^{-1}\va_i}_2^2}
	= \frac{(\va_i^\intercal(\mA^\intercal\mA)^{-1}\va_i)^2}{\va_i^\intercal(\mA^\intercal\mA)^{-1}\va_i}
	= \va_i^\intercal(\mA^\intercal\mA)^{-1}\va_i
\end{align}
Which completes the proof.
\end{proof}
\end{dropdown}

---

\begin{lemma}{Minimum Characterization}{min-char}
Let $\va_i$ be the $i^{th}$ row of $\mA$.
Then
$$
	\tau_i(\mA) = \min_{\vy\in\bbR^n,\,\mA^\intercal\vy=\va_i} \normof{\vy}_2^2
$$
\end{lemma}
\begin{dropdown}{_Proof_}
\begin{proof}
For simplicity, we assume that $\mA$ is both full-rank and tall-and-skinny.

Notice this minimization problem a minimum-norm underdetermined least-squares problem, with known solution $\vy=\mA(\mA^\intercal\mA)^{-1}\va_i$.
So, we know that
\[
	\min_{\vy\in\bbR^n,\,\mA^\intercal\vy=\va_i} \normof{\vy}_2^2
	= \normof{\mA(\mA^\intercal\mA)^{-1}\va_i}_2^2
	= \va_i^\intercal(\mA^\intercal\mA)^{-1}\va_i
	= \tau_i(\mA)
\]
where the second equality is shown in the start to the proof of \lemmaref{inner-product-char}, and the last equality is \lemmaref{inner-product-char}.
\end{proof}
\end{dropdown}

---

\begin{lemma}{Basic Properties of $\tau_i(\mA)$}{basic-props}
Let $\mA\in\bbR^{n \times d}$ full-rank with $n \geq d$.
Then,
1. Each leverage score has $\tau_i(\mA)\in[0,1]$
2. If $\mB\in\bbR^{d \times d}$ is full-rank, then $\tau_i(\mA\mB) = \tau_i(\mA)$
3. If $\mA^\intercal\mA=\mI$, then $\tau_i(\mA)=\normof{\va_i}_2^2$
4. If $\mU\in\bbR^{n \times d}$ with $\mU^\intercal\mU=\mI$ has the same columnspace as $\mA$, then $\tau_i(\mA)=\normof{\vu_i}_2^2$ where $\vu_i$ is the $i^{th}$ row of $\mU$.
5. The sum of leverages is $\sum_{i=1}^n \tau_i(\mA)=d$
\end{lemma}
\begin{dropdown}{_Proof_}
\begin{proof}

Point 1 follows directly from the Max Characterization (\definitionref{max-char}).

Point 2 also follows form the Max Characterization.
In particular, for any $\vx\in\bbR^{d}$ we can define $\vy\defeq\mB\vx$.
Since $\mB$ is invertible, the maximizing over all $\vx\in\bbR^d$ is equivalent to maximizing over all $\vy\in\bbR^d$:
\[
	\tau_i(\mA\mB)
	= \max_{\vx\in\bbR^d} \frac{[\mA\mB\vx]_i^2}{\normof{\mA\mB\vx}_2^2}
	= \max_{\vy\in\bbR^d} \frac{[\mA\vy]_i^2}{\normof{\mA\vy}_2^2}
	= \tau_i(\mA)
\]

Point 3 follows from the Inner Product Characterization:
\[
	\tau_i(\mA) = \va_i^\intercal(\mA^\intercal\mA)^{-1}\va_i = \va_i^\intercal\va_i = \normof{\va_i}_2^2
\]

Point 4 follows from Points 2 and 3.

Point 5 follows from Point 4:
Since every $\mA$ has an SVD $\mA=\mU\mat{\Sigma}\mV^\intercal$, we can find an orthogonal $\mU$ that satisfied Point 4.
Let $\vu_i$ denote the $i^{th}$ row of $\mU$ and let $\hat{\vu}_j$ denote the $j^{th}$ column of $\mU$.
Then,
\[
	\sum_{i=1}^n \tau_i(\mA)
	= \sum_{i=1}^n \normof{\vu_i}_2^2
	= \sum_{i=1}^n \sum_{j=1}^d \mU_{ij}^2
	= \sum_{j=1}^d \normof{\hat{\vu}_j}_2^2
	= d
\]
\end{proof}
\end{dropdown}

There are some intuitive implications from these bullet points:
- The leverage scores of $\mA$ depend _only_ on the range of $\mA$. Any other matrix with the same column space has the same leverage scores.
- If $\mA=\mQ\mR$ is the economic QR decomposition of $\mA$, and if $\vq_1,\ldots,\vq_n$ are the rows of $\mQ$, then $\tau_i(\mA)=\normof{\vq}_2^2$.

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

* \biblabel{alaoui2014fast}{Alaoui Mahoney (2015)} **Alaoui** and **Mahoney**. [Fast Randomized Kernel Methods with Statistical Guarantees.](https://www.stat.berkeley.edu/~mmahoney/pubs/elalaoui-nips15.pdf). _NIPS_ 2015.

* \biblabel{avron2011randomized}{Avron et al. (2017)} **Avron**, **Kapralov**, **Musco**, **Musco**, **Velingker**, and **Zandieh**. [Random Fourier Features for Kernel Ridge Regression: Approximation Bounds and Statistical Guarantees](https://arxiv.org/pdf/1804.09893.pdf). _ICML_ 2017.

* \biblabel{avron2019universal}{Avron et al. (2019)} **Avron**, **Kapralov**, **Musco**, **Musco**, **Velingker**, and **Zandieh**. [A Universal Sampling Method for Reconstructing Signals with Simple Fourier Transforms](https://arxiv.org/pdf/1812.08723.pdf). _STOC_ 2019.

* \biblabel{chen2019active}{Chen Price (2019)} **Chen** and **Price**. [Active Regression via Linear-Sample Sparsification](https://arxiv.org/pdf/1511.07263.pdf). _COLT_ 2019.

* \biblabel{avron2011input}{Cohen et al. (2017)} **Cohen**, **Kapralov**, **Musco**, and **Musco**. [Input Sparsity Time Low-Rank Approximation via Ridge Leverage Score Sampling](https://arxiv.org/pdf/1511.07263.pdf). _SODA_ 2017.

* \biblabel{cohen2015lp}{Cohen Peng (2019)} **Cohen** and **Peng**. [$\ell_p$ Row Sampling by Lewis Weights](https://arxiv.org/pdf/1412.0588.pdf). _STOC_ 2015.

\theoremscripts
