@def title = "RandNLA Proof Wiki"
@def tags = ["syntax", "code"]

\enabletheorems

# The Johnson-Lindenstrauss Lemma

The Johnsons-Lindenstrauss Lemma, or JL Lemma, says that we can embed $n$ points from $\bbR^d$ in a $O(\frac{\log(n)}{\eps^2})$ dimensional space while preserving all pairwise distances to relative error $\eps$.
We prove this in the simplest case, where we choose this embedding to be a random Gaussian matrix.
We start with a simple concentration for a single vector $\vx\in\bbR^{d}$:
\begin{lemma}{Gaussians Preserve Norms}{gaussian-norm}
Fix $\eps>0$, $\delta>0$, dimension $d$, and $k = \Omega(\frac1{\eps^2}\log(\frac1\delta))$.
Let $\mG\in\bbR^{k \times d}$ be a matrix of iid $\cN(0,1)$ entries.
Then with probability $1-\delta$, for any fixed $\vx\in\bbR^d$, we have
\[
	(1-\eps)\normof{\vx}_2^2 \leq \normof{\tsfrac{1}{\sqrt k}\mG\vx}_2^2 \leq (1+\eps)\normof{\vx}_2^2
\]
\end{lemma}
\begin{dropdown}{_Proof_}
\begin{proof}
Consider the $i^{th}$ entry of $\vw \defeq \mG\vx$.
Letting $g_1,\ldots,g_d$ denote the entries of the $i^{th}$ row of $\mG$, we have
\[
	w_i
	= \sum_{j=1}^d g_j x_j
	\sim \sum_{j=1}^d \cN(0,x_j^2)
	= \cN(0,\normof{\vx}_2^2)
\]
That is, $\vw$ is just a vector of iid $\cN(0,\normof{\vx}_2^2)$ random variables.
So,
\[
	\frac1k\normof{\vw}_2^2
	= \frac1k \sum_{i=1}^k w_i^2
	\sim \normof{\vx}_2^2 \cdot \frac1k\sum_{i=1}^k(\cN(0,1))^2
\]
Fortunately, the concentration of the average of squared Gaussians is well understood.
Equation 2.2.1 from \cite{wainright2019high} says that, with probability $1-\delta$, standard normal random variables $z_1,\ldots,z_k$ have
\[
	\abs{\frac{1}{k}\sum_{i=1}^k z_i^2 - 1} \leq \sqrt{\frac{8}{k}\ln(\frac2\delta)} \leq \eps
\]
where our value of $k = \Omega(\frac1{\eps^2}\log(\frac1\delta))$ guarantees this error term is at most $\eps$.
We conclude that
\[
	\abs{\normof{\vx}_2^2 - \normof{\tsfrac{1}{\sqrt k} \mG\vx}_2^2}
	= \abs{ \normof{\vx}_2^2 - \tsfrac{1}{k}\normof{\vw}_2^2 }
	= \normof{\vx}_2^2 \cdot \abs{ \frac{1}{k}\sum_{i=1}^k z_i^2 - 1 }
	\leq \eps \normof{\vx}_2^2
\]
Which completes the proof.
\end{proof}
\end{dropdown}

Given this simple concentration, we can state the Johnson-Lindenstrauss Lemma:

\begin{lemma}{Johnson-Lindenstrauss}{jl-lemma}
Let $\vx_1,\ldots,\vx_n\in\bbR^{d}$.
Fix $\eps\in(0,1)$, $\delta>0$, and $k = \Omega(\frac1{\eps^2}\log(\frac n\delta))$.
Let $\mPi\in\bbR^{k \times d}$ be a matrix of iid $\cN(0,\frac1k)$ entries.
Then, with probability $1-\delta$, for all pairs $i,j$ we have
\[
	(1-\eps) \normof{\vx_i - \vx_j}_2
	\leq \normof{\mPi\vx_i - \mPi\vx_j}_2
	\leq (1+\eps) \normof{\vx_i - \vx_j}_2
\]
\end{lemma}
\begin{dropdown}{_Proof_}
\begin{proof}
For all pairs $i,j$, consider $\vv_{i,j} \defeq \vx_i - \vx_j$, and union bound \lemmaref{gaussian-norm} to find a matrix $\mG$ which preserves the norms of all $\vv_{i,j}$.
This involves union bounding over $\binom{n}{2} = O(n^2)$ vectors, which causes $\log(\frac n\delta)$ to appear in the requirement on $k$.
Since $\mPi = \frac1{\sqrt k} \mG$, we have
\[
	\sqrt{1-\eps} \normof{\vx_i - \vx_j}_2
	\leq \normof{\mPi\vx_i - \mPi\vx_j}_2
	\leq \sqrt{1+\eps} \normof{\vx_i - \vx_j}_2
\]
We then note that $\sqrt{1-\eps} > 1-\eps$ and $\sqrt{1+\eps} < 1+\eps$ ([see this on Desmos](https://www.desmos.com/calculator/heb771c0tf)), which completes the proof.
\end{proof}
\end{dropdown}

# See Also

The proofs above are ubiquitous, but can be found for example in \cite{wainright2019high} and \cite{musco18lecture}.

Here's some important papers in the area of JL:
- \cite{johnson1984extensions} is the original paper of Johnson and Lindenstrauss.
- \cite{musco18lecture} has a nice short proof, which this page basically copies.
- \cite{larsen2017optimality} shows that no other mapping, even possibly nonlinear mapping, can outperform the rate of $k = \Omega(\frac{1}{\eps^2} \log(n))$.
- \cite{ailon2009fast} shows the Fast-JL Embedding, which samples columns from Fourier or Hadamard matrices, allowing for faster computation of the matrix-vector product $\mPi\vx$.
- \cite{kane2014sparser} gives the sparsest known embedding matrix $\mPi$.
- \cite{nelson2013sparsity} lower bounds the neccessary sparsity.
- _Let me know if anything is missing_

# Bibliography

* \biblabel{ailon2009fast}{Ailon Chazelle (2009)} **Ailon** and **Chazelle**. [The fast Johnson--Lindenstrauss transform and approximate nearest neighbors](https://www.cs.princeton.edu/~chazelle/pubs/FJLT-sicomp09.pdf). _SICOMP 2009_.

* \biblabel{johnson1984extensions}{Johnson Lindenstrauss (1984)} **Johnson** and **Lindenstrauss**. [Extensions of Lipschitz mappings into a Hilbert space](http://stanford.edu/class/cs114/readings/JL-Johnson.pdf). _Contemporary Mathematics 1984_.

* \biblabel{kane2014sparser}{Kane Nelson (2014)} **Kane** and **Nelson**. [Sparser Johnson-Lindenstrauss Transforms](https://arxiv.org/pdf/1012.1577.pdf). _JACM 2014_.

* \biblabel{larsen2017optimality}{Larsen Nelson (2017)} **Larsen** and **Nelson**. [Optimality of the Johnson-Lindenstrauss Lemma](https://arxiv.org/pdf/1609.02094.pdf). _FOCS 2017_.

* \biblabel{musco18lecture}{Musco (2018)} **Musco**. [Lecture 10: Dimensionality Reduction and the Johnson-Lindenstrauss Lemma](https://www.cs.princeton.edu/courses/archive/fall18/cos521/Lectures/lec10.pdf). _Lecture Notes, 2018_.

* \biblabel{nelson2013sparsity}{Nelson Nguyễn (2013)} **Nelson** and **Nguyễn**. [Sparsity Lower Bounds for Dimensionality Reducing Maps](https://arxiv.org/pdf/1211.0995.pdf). _STOC 2013_.

* \biblabel{wainright2019high}{Wainwright (2015)} **Wainwright**. [Draft of _High-dimensional statistics: A Non-Asymptotic Viewpoint_, Chapter 2](https://www.stat.berkeley.edu/~mjwain/stat210b/Chap2_TailBounds_Jan22_2015.pdf). _Draft of publication at Cambridge University Press, 2015_.

\theoremscripts
