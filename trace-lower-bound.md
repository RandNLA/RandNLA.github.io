@def title = "Trace Estimation Lower Bounds (Hidden Wishart Method)"

\enabletheorems

# Introduction

Consider the problem of estimating the trace of a PSD matrix $\mA$ from a small number of matrix-vector products.
There's a variety of algorithms which achieve this goal using merely $\cO(1/\eps)$ matrix-vector products.
Here, we will show a succint proof of the complementary lower bound: that $\Omega(1/\eps)$ matrix-vector products is necessary in the worst case.

\begin{theorem}{Trace Estimation Lower Bound}{trace-lb}
    Any algorithm that accesses a PSD matrix $\mA$ via $k$ (possibly adaptive) matrix-vector products and outputs an estimate $\tilde{t}$ of $\tr(\mA)$ such that $\sqrt{\E[|\tilde{t} - \tr(\mA)|^2]} \leq \eps\tr(\mA)$ must use at least $k \geq \frac{1}{2\sqrt{2}\eps}$ matrix-vector products.
\end{theorem}

To prove this lower bound, we will need two ingredients: the hidden Wishart theorem, and the conditional expectation.

## The Hidden Wishart Theorem

We will rely on the following remarkable result, whose proof we omit.
This is the core theorem that enables this entire lower bound technique.

\begin{theorem}{Hidden Wishart Theorem}{hidden-wishart}
    Let $\mG \in \bbR^{n \times n}$ be a random matrix with iid $\cN(0,1)$ entries, and let $\mA = \mG^\intercal\mG$.
    Suppose an algorithm computes $k$ (possibly adaptive) matrix-vector products with $\mA$, denoted $\vy_1 = \mA \vx_1, \ldots, \vy_k = \mA \vx_k$.

    Then, there exists a matrix $\mDelta \in \bbR^{n \times n}$ and orthogonal matrix $\mV \in \bbR^{n \times n}$, each constructed deterministically from the queries $\{(\vx_i, \vy_i)\}_{i=1}^k$, such that
    \[
        \mV^\intercal\mA\mV = \mDelta + \begin{bmatrix} \mat{0} & \mat{0} \\ \mat{0} & \tilde{\mA} \end{bmatrix},
    \]
    where $\tilde{\mA}=\tilde{\mG}^\intercal\tilde{\mG} \in \bbR^{(n-k) \times (n-k)}$ and $\tilde{\mG} \in \bbR^{(n-k) \times (n-k)}$ has iid $\cN(0,1)$ entries independent of $\{(\vx_i, \vy_i)\}_{i=1}^k$.
\end{theorem}

A succinct proof of \theoremref{hidden-wishart} can be found in Appendix B.1 of \cite{amsel26}.

The matrices $\mA$ and $\tilde{\mA}$ follow the _Wishart_ distribution, hence the name of the theorem.
The result shows that after $k$ matrix-vector products, there remains a large random component $\tilde{\mA}\in\bbR^{(n-k) \times (n-k)}$ of the matrix $\mA\in\bbR^{n \times n}$ which is completely independent of the algorithm's matrix-vector queries.
This is the titular "hidden Wishart".

Notice that this theorem has robbed the matrix-vector algorithm of any and all agency -- no matter how the method chooses its (possibly adaptive) queries, there is always a large random component of the matrix which it has no information about.
Since we can characterize a large part of $\mA$ that the algorithm has no information about at all, we can use simple statistical tools to prove \theoremref{trace-lb}.

## A simple statistical observation

In \theoremref{trace-lb}, we care about minimizing the mean squared error of our algorithms given some data about $\mA$ (namely, a sequence of matrix-vector products).
Classical statistics tells us that the best possible algorithm for this goal is the conditional expectation:

\begin{lemma}{Conditional Expectation Minimizes MSE}{mmse}
    Let $X$ and $Y$ be (possibly dependent) random variables.
    Suppose an algorithm observes $Y$ and outputs an estimate $\tilde{X}$ of $X$ based on $Y$.
    Then, the error of $\tilde{X}$ is lower bounded as
    \[
        \E\big[|\tilde{X} - X|^2\big] \geq \E\big[\Var[X\,|\,Y]\big],
    \]
    and this lower bound is achieved by the conditional expectation $\hat{X} \defeq \E[X\,|\,Y]$.
\end{lemma}
\begin{dropdown}{_Proof_}
\begin{proof}
    By the tower rule,
    \[
        \E\big[|\tilde{X} - X|^2\big] = \E\left[\E\big[|\tilde{X} - X|^2 ~|~ Y\big]\right].
    \]
    For any fixed $Y=y$, the inner expectation $\E\big[|\tilde{X} - X|^2 ~|~ Y=y\big]$ is minimized by choosing $\tilde{X} = \E[X\,|\,Y=y]$.
    Depending on who you ask, this is either the definition of conditional expectation or a basic fact about it.
    Either way, we have
    \[
        \E\big[|\tilde{X} - X|^2 ~|~ Y\big] \geq \E\left[\big| X - \E[X|Y] \big|^2 ~|~ Y\right] = \Var[X\,|\,Y].
    \]
    Taking the expectation over $Y$ finished the proof.
\end{proof}
\end{dropdown}


## Proof of Trace Estimation Lower Bound

From this result, we can now prove the theorem.
It will be helpful to keep in mind that if $Z \sim \chi^2_d$ is a chi-squared random variable with $d$ degrees of freedom, then $\E[Z] = d$ and $\Var[Z] = 2d$.

\begin{dropdown}{_Proof of \theoremref{trace-lb}_}
\begin{proof}
    Let $\mA$ be defined as in \theoremref{hidden-wishart}.
    Note that $\E[\tr(\mA)] = \E[\tr(\mG^\intercal\mG)] = \E[\normof{\mG}_{\rm F}^2] = n^2$.
    By \theoremref{hidden-wishart}, after $k$ matrix-vector products, there exists a decomposition
    \[
        \mV^\intercal\mA\mV = \mDelta + \begin{bmatrix}\mat{0} & \mat{0} \\ \mat{0} & \tilde{\mA}\end{bmatrix},
    \]
    where $\tilde{\mA}\in\bbR^{(n-k) \times (n-k)}$ is independent of the algorithm's queries.
    By \lemmaref{mmse}, the lowest possible error any algorithm can achieve is $\E[\Var[\tr(\mA) \mid \mDelta,\mV]]$.
    We note that
    \[
        \tr(\mA) = \tr(\mDelta) + \tr(\tilde{\mA}).
    \]
    Since $\tilde{\mA}$ is independent of $\mDelta$ and $\mV$, and since $\tr(\tilde{\mA})=\normof{\tilde{\mG}}_{\rm F}^2$ has a $\chi^2$ distribution with $(n-k)$ degrees of freedom, we have
    \[
        \E[|\tilde{t} - \tr(\mA)|^2] \geq \Var[\tr(\mA) \mid \mDelta,\mV] = \Var[\tr(\tilde{\mA})] = 2(n-k)^2.
    \]
    So, any estimator that achieves root mean squared error at most $\eps \tr(\mA)$ must satisfy
    \[
        \sqrt{2}(n-k) \leq \sqrt{\E[|\tilde{t} - \tr(\mA)|^2]} \leq \eps \E[\tr(\mA)] = \eps n^2.
    \]
    Rearranging this inequality yields
    \[
        k \geq n - \frac{\eps n^2}{\sqrt{2}}.
    \]
    Maximizing the right-hand side over $n$ gives the desired lower bound of $k \geq \frac{1}{2\sqrt{2}\eps}$.
\end{proof}
\end{dropdown}

# See Also

There have been many papers that use the Hidden Wishart Theorem to prove matrix-vector complexity lower bounds.
The proof here is a special case of an analysis in \citep{meyer2023}.

Other relevant papers to this article include:

- \cite{braverman2020} Introduces the hidden Wishart method to prove lower bounds for linear regression and eigenvalue estimation.
- \cite{meyer2021} Has the first optimal $\Omega(1/\eps)$ lower bounds for trace estimation, but uses more complex methods.
- \cite{amsel26} Gives a succinct proof of the hidden Wishart theorem.
- \cite{jiang21} Sharpens the analysis above to get a nearly tight dependence on failure probability.
- \cite{amsel25} Uses the hidden Wishart method to prove lower bounds for learning the diagonal of a matrix.
- _Let me know if anything is missing. I'm sure many papers are missing._

# References

* \biblabel{braverman2020}{Braverman et al. (2020)} **Braverman**, **Hazan**, **Simchowitz**, and **Woodworth**. [The gradient complexity of linear regression](https://arxiv.org/pdf/1911.02212v3). _COLT_ 2020.
* \biblabel{jiang21}{Jiang et al. (2021)} **Jiang**, **Pham**, **Woodruff**, and **Zhang**. [Optimal Sketching for Trace Estimation](https://arxiv.org/pdf/2111.00664v1). _NeurIPS_ 2021.
* \biblabel{meyer2021}{Meyer et al. (2021)} **Meyer**, **Musco**, **Musco**, and **Woodruff**. [Hutch++: Optimal stochastic trace estimation](https://arxiv.org/pdf/2010.09649v5). _SOSA_ 2021.
* \biblabel{meyer2023}{Meyer Avron (2023)} **Meyer** and **Avron**. [Hutchinson's Estimator is Bad at Kronecker-Trace-Estimation](https://arxiv.org/pdf/2309.04952v2). _preprint_ 2023.
* \biblabel{amsel25}{Amsel et al. (2025)} **Amsel**, **Avi**, **Chen**, **Keles**, **Hegde**, **Musco**, **Musco**, and **Persson**. [Query Efficient Structured Matrix Learning](https://www.arxiv.org/pdf/2507.19290v1). _preprint_ 2025.
* \biblabel{amsel26}{Amsel et al. (2026)} **Amsel**, **Chen**, **Keles**, **Halikias**, **Musco**, and **Musco**. [Fixed-sparsity matrix approximation from matrix-vector products](https://arxiv.org/pdf/2402.09379v3). _SIMAX_ 2026.

<!-- The best constants I know for the upper bound come from using Tropp Webber for Nystrom++ and gets about 7.84/eps matvecs -->

\theoremscripts
