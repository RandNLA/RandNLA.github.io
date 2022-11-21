@def title = "RandNLA Proof Wiki"
@def tags = ["syntax", "code"]

\enabletheorems

# Subspace Embedding via $\eps$-Nets

On this page we prove that a Johnson-Lindenstrauss Transform is an \highlight{oblivious subspace embedding}.
That is, the JL matrix is a map from $\bbR^n$ to $\bbR^{O(d)}$ which preserves the norms of all vectors in any fixed $d$-dimensional subspace of $\bbR^n$.
This is called _oblivious_ because the JL map doesn't know anything about the subspace beyond how big it is (i.e. it's dimension).

_Prerequisite: [Gaussian Johnson-Lindenstrauss Lemma](/johnson-lindenstrauss/)_

\begin{theorem}{Subspace Embedding}{subspace-embed}
Fix $\eps\in(0,1)$ and $\delta\in(0,1)$.
Let $\cV$ be a $d$-dimensional linear subspace of $\bbR^n$.
Let $\mPi\in\bbR^{k \times n}$ be a JL sketch for $k = \Omega(\frac{d + \log(1/\delta)}{\eps^2})$.
Then with probability $1-\delta$, for all $\vx\in\cV$ we have
\[
	(1-\eps)\normof{\vx}_2 \leq \normof{\mPi\vx}_2 \leq (1+\eps) \normof{\vx}_2
\]
\end{theorem}

We will prove this via an \highlight{$\eps$-Net Argument}, which works in two stages.
First, we build a "Net" of $O( (\frac1\eps)^d)$ many points and union-bound the JL guarantee over all of these fixed points.
Second, we make a \highlight{rounding argument}, which says that since linear maps are smooth, the JL guarantee on the net implies the JL guarantee on all points in $\cV$.

<!--
Intuitively, we can think of the net as discretizing some important set (e.g. the unit ball), and the rounding argument says that if the discretization is fine enough then the worst error on the whole set equals the worst error on the discrete points.
-->

We show two different rounding arguments.
The first is simpler and just uses the triangle inequality, but incurs a suboptimal $k=\Omega(\frac{d\log(1/\eps) + \log(1/\delta)}{\eps^2})$ dependence.
The second is more careful and leverages the relationship between the $\ell_2$ norm and the inner product, which achieves the rate in \theoremref{subspace-embed}.

# Building a Net

Our first step in the proof is to build an fine net.
Note that, by linearity, the guarantee from \theoremref{subspace-embed} is equivalent to saying that all unit vectors in $\cV$ have their norms preserved:
\[
	\normof{\mPi\vx}_2 \in (1\pm \eps)
	\hspace{1cm}
	\forall \vx\in\cV ~ \text{ s.t. } \normof{\vx}_2=1
\]
That is, we only need to verify that $\normof{\mPi\vx}_2$ is approximately accurate for the unit sphere in $\cV$.
So, we will only build a net over the unit sphere in $\cV$:

\begin{theorem}{$\ell_2$ Net Size}{net-size}
Fix $\eps\in(0,2)$, and dimension $d\geq1$.
Then there exists a set $\cN_\eps$ with $\abs{\cN_\eps} \leq (\frac 6\eps)^d$ such that, for all $\vx$ with $\normof{\vx}_2=1$ there exists some $\vy\in\cN_\eps$ such that $\normof{\vx-\vy}_2\leq\eps$.
\end{theorem}
\begin{dropdown}{_Proof via Volume Argument_}
\begin{proof}
Consider a greedy algorithm that constructs $\cN_\eps$:
- Start with $\cN_\eps=\{\}$
- While there exists any $\vx$ with $\normof{\vx}_2=1$ such that $\normof{\vx-\vy}_2>\eps$ for all $\vy\in\cN$, add $\vx$ to $\cN_\eps$

Then $\cN_\eps$ clearly satisfies the correctness property -- all $\vx$ on the unit ball must be $\eps$-close to some $\vy\in\cN_\eps$.
We just have to bound $\abs{\cN_\eps}$.
Note that for all $\vy_i,\vy_j\in\cN_\eps$, we have $\normof{\vy_i-\vy_j}_2>\eps$, or else one of those $\vy$ vectors would not have been added by the greedy algorithm.

Note that a ball of radius $r$ in $\bbR^d$ has volume $cr^d$, where $c = \frac{\pi^{d/2}}{\Gamma(\frac{d}{2}+1)}$ does not depend on $r$ ([link](https://en.wikipedia.org/wiki/Volume_of_an_n-ball#Formulas)).

We now make the \highlight{Volume Argument}.
We consider placing balls of radius $\frac\eps2$ on all $\vy_i\in\cN_\eps$:
\[
	\cB(\vy_1,\tsfrac{\eps}{2}), ~ \ldots, ~ \cB(\vy_{\abs{\cN_\eps}},\tsfrac{\eps}{2})
\]
Since these balls do not intersect, and since the volume of each ball is $c\, (\frac{\eps}{2})^d$, these balls have total volume equal to $\abs{\cN_\eps} \, c \, (\frac{\eps}{2})^d$.
Next, since the _union_ of these disjoint balls all fits within a ball of radius $1+\frac{\eps}{2} \leq 3$, we must have
\[
	\abs{\cN_\eps} \, c \, (\tsfrac{\eps}{2})^d \leq c \, 3^d
\]
Rearranging this, we get $\abs{\cN_\eps} \leq (\frac{6}{\eps})^d$.
\end{proof}
\end{dropdown}
\begin{dropdown}{Comments on the Proof}
- \theoremref{net-size} does not care what $d$-dimensional space is considered, be it $\bbR^d$ or a $d$-dimensional linear subspace of $\bbR^n$.
We will use it to make nets over the unit ball of $\cV$.
- The _Volume Argument_ used in the proof does not really use the fact that $\normof{\vx}_2=1$, and instead actually bounds the size of net over the whole ball of $\normof{\vx}_2\leq1$.
So, if you need an argument that uses a net over the whole interior of the ball, the same bound of $\abs{\cN_\eps} \leq (\frac6\eps)^d$ is correct.
- The same proof when constrained to a more typical $\eps\in(0,1)$ achieves a different constant $\abs{\cN_\eps}\leq (\frac4\eps)^d$, which may appear in other works.
\end{dropdown}

# Net Expansion of a Vector

<!--
We now show a boilerplate method for proving \theoremref{subspace-embed}, using just the triangle inequality.
We union bound the JL guarantee over all vectors $\vy\in\cN_\eps$, so we know that $\normof{\mPi\vy}_2\approx\normof{\vy}_2$ for everything in the net.
To prove \theoremref{subspace-embed}, as we discussed, we have to prove by the smoothness of the matrix  $\normof{\mPi\vx}_2\approx 1$ for all $\normof{\vx}_2=1$.
To achieve this, we expand any such $\vx$ as a converging sum of items in the net:
-->

To prove \theoremref{subspace-embed}, we first have to represent an arbitrary $\vx$ on the unit ball in terms of points on the net:

\begin{lemma}{Net Expansion of a Vector}{net-expansion}
Let $\normof{\vx}_2=1$ and let $\cN_\eps$ be an $\eps$-Net for the unit ball.
Then, there exists a sequence $\vy_0,\ldots,\vy_n,\ldots\in\cN_\eps$ such that $\vx=\sum_{i=0}^\infty \alpha_i \vy_i$ where $0 \leq \alpha_i \leq \eps^i$ and $\alpha_0=1$.
\end{lemma}
\begin{dropdown}{_Proof_}
\begin{proof}
By construction of the net, we know there exists some $\vy_0\in\cN_\eps$ such that $\normof{\vx-\vy_0}\leq\eps$.
That is, the residual $\vr_0 \defeq \vx-\vy_0$ has norm $c_1\defeq\normof{\vr_0}\leq\eps$.
Then, again by the net, we know that some $\vy_1\in\cN_\eps$ has $\normof{\frac{\vr_0}{c_1} - \vy_1} \leq \eps$.
That is, the residual $\vr_1 \defeq \frac{\vr_0}{c_1} - \vy_1$ has norm $c_2\defeq\normof{\vr_1}\leq\eps$.
Repeating this process, we get
\begin{align}
	\vx
	&= \vy_0 + \vr_0 \\
	&= \vy_0 + c_1(\vy_1 + \vr_2) \\
	&= \vy_0 + c_1(\vy_1 + c_2(\vy_2 + \ldots)) \\
	&= \vy_0 + c_1\vy_1 + c_1c_2\vy_2 + c_1c_2c_3\vy_3 + \ldots
\end{align}
Since each $c_i\in[0,\eps]$, we get that $\alpha_i \defeq c_1c_2\ldots c_i \in[0,\eps^i]$.
\end{proof}
\end{dropdown}

<!--
We are now ready to complete net the argument, starting with the simpler method with slightly worse rates.
-->

# Rounding via Triangle Inequality

Here, we present a proof of \theoremref{subspace-embed}, but which uses sample complexity $O(\frac{d \log(1/\eps) + \log(1/\delta)}{\eps^2})$ instead of the tighter rate $O(\frac{d + \log(1/\delta)}{\eps^2})$ promised in theorem statement.
This proof, however, is very simple and just uses \lemmaref{net-expansion} and the triangle inequality:
\begin{dropdown}{_Proof of \theoremref{subspace-embed} for $k=\Omega(\frac{d\log(1/\eps) + \log(1/\delta)}{\eps^2})$_}
\begin{proof}
Let $\eps_0 \defeq \frac{\eps}{4}$.
Let $\cN_{\eps_0}$ be an $\eps_0$-Net for the unit ball, so that $\abs{\cN_{\eps_0}}\leq (\frac{24}\eps)^d$ and so that union bounding JL over all $\vy\in\cN_{\eps_0}$ requires sketching dimension $k=\Omega(\frac{\log(\abs{\cN_{\eps_0}}/\delta)}{\eps^2}) = \Omega(\frac{d \log(1/\eps) + \log(1/\delta)}{\eps^2})$.
Namely, we have $\normof{\mPi\vy_i}\in(1\pm\eps_0)$ since $\normof{\vy}_2=1$ for all $\vy\in\cN_{\eps_0}$.

Let $\vx$ be any vector with $\normof{\vx}_2=1$, and let $\vx = \sum_{i=0}^\infty \alpha_i \vy_i$ be its Net Expansion.
Then, we have
\[
	\normof{\mPi\vx}_2
	\leq \sum_{i=0}^\infty \alpha_i\normof{\mPi\vy_i}
	\leq (1+\eps_0)\sum_{i=0}^\infty \eps_0^i
	= \frac{1+\eps_0}{1-\eps_0}
\]
since $\frac{1+\eps_0}{1-\eps_0} \leq 1+4\eps_0 = 1+\eps$ for $\eps\in[0,1]$, we get $\normof{\mPi\vx}_2\leq1+\eps$.
Similarly, by the reverse triangle inequality $\normof{\va+\vb}\geq\normof{\va}-\normof{\vb}$,
\begin{align}
	\normof{\mPi\vx}
	&\geq \normof{\mPi\vy_0} - \normoflr{\sum_{i=1}^\infty \alpha_i\mPi\vy_i} \\
	&\geq \normof{\mPi\vy_0} - \sum_{i=1}^\infty \alpha_i \normof{\mPi\vy_i} \\
	&\geq (1-\eps_0) - (1+\eps_0) \sum_{i=1}^\infty \eps_0^i \\
	&= (1-\eps_0) - (1+\eps_0) \frac{\eps_0}{1-\eps_0} \\
\end{align}
So we get $\normof{\mPi\vx} \geq (1-\eps_0) - \eps_0 \frac{1+\eps_0}{1-\eps_0} \geq 1-3\eps_0 \geq 1-\eps$.
That is, we overall find
\[
	\normof{\mPi\vx}_2 \in (1\pm \eps)
	\hspace{1cm}
	\forall \vx\in\cV ~ \text{ s.t. } \normof{\vx}_2=1
\]
Which completes the proof.
\end{proof}
\end{dropdown}

# Rounding via Inner Products

We now present a sharper analysis that achieves the rate of $k=\Omega(\frac{d+\log(1/\delta)}{\eps^2})$:

\begin{dropdown}{_Intuition for the Proof_}
We do this by decreasing the precision of the net from an $\eps$-Net to a $\frac12$-Net, which therefore needs a new tighter rounding argument.
To understand how this works, we ask why the proof via triangle inequality has to use a net with precision $O(\eps)$.
Suppose the JL matrix preserved the norm of all $\vy\in\cN_{\eps}$ perfectly, then that proof bounds
\[
	\normof{\vx}\leq\sum_{i=0}^\infty \alpha_i \normof{\vy_i} \leq \sum_{i=0}^\infty \eps^i = 1+\eps
\]
This proof, even for a perfectly accurate JL matrix, still overestimates the norm of $\vx$ because the triangle inequality is loosing vital information.
Specifically, note that $\normof{\va+\vb}_2 = \normof{\va}_2 + \normof{\vb}_2$ only if $\va$ and $\vb$ are perfectly orthogonal.
If we somehow had a net $\cN_\eps$ such that $\vy_1,\ldots,\vy_\infty$ were orthogonal, then the triangle inequality proof would be tight.

However, that's trivially not the case here.
For instance, by the pigeonhole principle we guarantee that $\vy_1,\ldots,\vy_\infty$ has infinitely many repeated vectors, since $\abs{\cN_\eps}<\infty$, and so those repeated vectors are deeply non-orthogonal.

So, we need to find a new way to preserve $\normof{\vx}$ in terms of $\vy_1,\ldots,\vy_\infty$, and that approach follows by examining the unique properties of the $\ell_2$ norm, namely that
\[
	\normof{\va-\vb}_2^2 = \normof{\va}_2^2 + \normof{\vb}_2^2 - 2\va^\intercal\vb
\]
Or, equivalently,
\[
	\va^\intercal\vb = \frac12 \left(\normof{\va}_2^2 + \normof{\vb}_2^2 - \normof{\va-\vb}_2^2\right)
\]
We then can preserve all three terms on the right to relative error by union bounding JL over $\va$, $\vb$, and $\va-\vb$, and 
So, we can expand $\normof{\vx}_2^2 = (\sum_i \alpha_i \vy_i)^\intercal(\sum_i \alpha_i \vy_i)$ as a large sum of inner products, preserve all the corresponding norms by JL, and recover a relative error guarantee with a coarser net and $(1+\eps)$-accurate JL.
\end{dropdown}
\begin{dropdown}{_Proof of \theoremref{subspace-embed}_}
Let $\eps_0 = \frac{\eps}{24}$.
Let $\cN_2$ be a $\frac12$-Net for the unit ball, so that $\abs{\cN_2} \leq 12^d$.
We union bound JL over all pairs $\vy,\vy'\in\cN_2$, so that both $\normof{\mPi\vy}_2\in(1\pm\eps_0)$ and $\normof{\mPi(\vy-\vy')}_2\in(1\pm\eps_0)\normof{\vy-\vy'}_2$ hold for all $\vy,\vy'\in\cN_2$.
This requires sketching dimension $k = \Omega(\frac{\log(\abs{\cN_2}/\delta)}{\eps_0^2}) = \Omega(\frac{d + \log(1/\delta)}{\eps^2})$.

Let $\vx$ be any vector with $\normof{\vx}_2=1$, and let $\vx = \sum_{i=0}^\infty \alpha_i \vy_i$ be its Net Expansion.
Then, we have
\begin{align*}
	\normof{\mPi\vx}_2^2
	&= (\textstyle{\sum_{i=0}^\infty} \alpha_i \mPi\vy_i)^\intercal (\textstyle{\sum_{j=0}^\infty} \alpha_j \mPi\vy_j) \\
	&= \sum_{i=0}^\infty \sum_{j=0}^\infty \alpha_i \alpha_j (\mPi\vy_i)^\intercal(\mPi\vy_j)
\end{align*}
We then bound the accuracy of this inner product, using the fact that $(1+\eps_0)^2 \leq 1+3\eps_0$ for $\eps_0 \leq 1$:
\begin{align*}
	(\mPi\vy_i)^\intercal(\mPi\vy_j)
	&= \frac12 \left( \normof{\mPi\vy_i}_2^2 + \normof{\mPi\vy_j}_2^2 - \normof{\mPi(\vy_i-\vy_j)}_2^2 \right) \\
	&\leq \frac12 \left( \normof{\vy_i}_2^2 + \normof{\vy_j}_2^2 - \normof{\vy_i-\vy_j}_2^2 \right) + \frac{3\eps_0}2 \left( \normof{\vy_i}_2^2 + \normof{\vy_j}_2^2 + \normof{\vy_i-\vy_j}_2^2 \right) \\
	&\leq \frac12 \left( \normof{\vy_i}_2^2 + \normof{\vy_j}_2^2 - \normof{\vy_i-\vy_j}_2^2 \right) + \frac{3\eps_0}2 \left( 1 + 1 + 2 \right) \\
	&= \vy_i^\intercal\vy_j + 6\eps_0
\end{align*}
With the matching lower bound following from $(1-\eps_0)^2 \geq 1-3\eps_0$.
So the subspace embedding error grows as
\begin{align*}
	\normof{\mPi\vx}_2^2
	&= \sum_{i=0}^\infty \sum_{j=0}^\infty \alpha_i \alpha_j (\mPi\vy_i)^\intercal(\mPi\vy_j) \\
	&\leq \sum_{i=0}^\infty \sum_{j=0}^\infty \alpha_i \alpha_j (\vy_i^\intercal\vy_j + 6\eps_0) \\
	&\leq \normof{\vx}_2^2 + \sum_{i=0}^\infty \sum_{j=0}^\infty \alpha_i \alpha_j 6\eps_0 \\
	&\leq 1 + 6\eps_0 \sum_{i=0}^\infty \sum_{j=0}^\infty 2^{-i} 2^{-j} \\
	&= 1 + 24\eps_0
\end{align*}
And the lower bound similarly is $\normof{\mPi\vx}_2^2 \geq 1 - 24\eps_0$.
So, we have $\normof{\mPi\vx}_2 \leq \sqrt{1+24\eps_0} \leq 1+24\eps_0 = 1+\eps$ and $\normof{\mPi\vx}_2 \geq \sqrt{1-24\eps_0} \leq 1-24\eps_0 = 1-\eps$ ([see this on Desmos](https://www.desmos.com/calculator/heb771c0tf)), which completes the proof.
\end{dropdown}

# See Also

The proofs above are ubiquitous, for example the coarser argument is in \cite{musco18lecture}.

Here's some important papers in world of oblivious subspace embeddings:
- \cite{johnson1984extensions} is the original paper of Johnson and Lindenstrauss.
- \cite{musco18lecture} has a nice short proof, which this page basically copies.
- _I'm sure **a lot** of great references are missing_
- _Let me know if anything is missing_

# Bibliography

* \biblabel{johnson1984extensions}{Johnson Lindenstrauss (1984)} **Johnson** and **Lindenstrauss**. [Extensions of Lipschitz mappings into a Hilbert space](http://stanford.edu/class/cs114/readings/JL-Johnson.pdf). _Contemporary Mathematics 1984_.

* \biblabel{musco18lecture}{Musco (2018)} **Musco**. [Lecture 11: Approximate regression, $\eps$-nets, and faster JL embeddings](https://www.cs.princeton.edu/courses/archive/fall18/cos521/Lectures/lec11.pdf). _Lecture Notes, 2018_.

\theoremscripts
