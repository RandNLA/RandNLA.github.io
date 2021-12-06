@def title = "RandNLA Proof Wiki"
@def tags = ["syntax", "code"]

\newcounter{NumAlgorithms}

# Contributing to the Wiki

We are more than happy to have people submit more proofs to the wiki, and this page explains the ideology behind the proofs, their style, and the technical requirements in submitting a proofs.
First, in ideology, our goal is to have a _minimal_ list of _essential_ and _simplified_ proofs in RandNLA.
So, we are especially interested in:
1. Proofs with a fundamentally different structure from what is already presented
2. Simplifications of proofs already on this site

If you have such a proof in mind, then please submit it!

There are three basic ways you can contribute:
1. Email us with a link/outline of a proof that should be on this page. We will respond and write this up when time permits, but this may take a while.
1. Email us a `.tex` document of the proof using our [LaTeX Template](/assets/RandNLA Wiki LaTeX Template.zip).
1. Submit a pull request to our git with your changes

The first approach may take a long time, and the fourth approach works best for small changes.
Since this website is written in Julia using Franklin.jl, and has some small custom scripts, submitting a pull request with a large change (e.g. a full new proof) may be more effort than is needed.
So we recommend the second and third approaches for large changes, especially new pages, as this lets us verify that all the flags are set correctly on the page.

In particular our LaTeX stylesheet is designed to be easy to merge into this website.
We list the relevant notation and definitions of the stylesheet below.
Once your proof is ready, submit it to Raphael Meyer at [ram900@nyu.edu](ram900@nyu.edu).

# Stylesheet API

For math commands, use:

| Meaning | Syntax | Appearance |
| --- | --- | --- |
| Matrices | `\mA, \mB, \mC, ..., \mZ` | $\mA$, $\mB$, $\mC$, ..., $\mZ$ |
| Vectors | `\va, \vb, \vc, ..., \vz` | $\va$, $\vb$, $\vc$, ..., $\vz$ |
| Scalars | `a, b, c, ..., z` | $a$, $b$, $c$, ..., $z$ |
| Black-board Bold | `\bbA, \bbB, \bbC, ..., \bbZ` \nbsp\nbsp\nbsp\nbsp\nbsp\nbsp | $\bbA$, $\bbB$, $\bbC$, ..., $\bbZ$ |
| Calligraphic Letters | `\cA, \cB, \cC, ..., \cZ` | $\cA$, $\cB$, $\cC$, ..., $\cZ$ |
| Probability | `\E[X], \Var[X], \Pr[A]` | $\E[X]$, $\Var[X]$, $\Pr[A]$ |
| Math Definition | `\defeq` | $\defeq$ |
| Epsilon | `\eps` | $\eps$ |
| Absolute Value | `\abs{...}` | $\abs{...}$ |
| Norm (without subscript) \nbsp\nbsp\nbsp\nbsp\nbsp\nbsp | `\normof{\mA}` | $\normof{\mA}$ |
| Norm (with subscript) | `\onormof{\mA}{2}` | $\onormof{\mA}{2}$ |
| Transpose | `\mA^\intercal\vx` | $\mA^\intercal\vx$ |
| Indicator Variable | `\indicate{i=j}` | $\indicate{i=j}$ |

If you want to use a nonstandard matrix or vector, such as a Greek letter matrix or vector, then you can type `\mat{\Lambda}` or `\vecalt{\alpha}` to get $\mat{\Lambda}$ or $\vecalt{\alpha}$.

For definitions, lemmas, and other standard `amsthm` environments, LaTeX users proceed as usual:
```
\begin{definition}[Some Name] %% Could be lemma, theorem, corollary, etc.
\label{named-label}
A definition is here. This is a good place to define things.
\end{definition}

\begin{proof}
This is a proof. Prove things here. This should probably follow a theorem, lemma, or corollary.
This relates to \definitionref{named-label}.
\end{proof}

\begin{algorithm}[Algo Name]
This is actually just like the lemma and theorem environments.
It is not tied to any of the standard latex pseudocode libraries.
This is intentional.
\end{algorithm}
```
Which yields the following:
\begin{definition}{Some Name}{named-label}
A definition is here. This is a good place to define things.
\end{definition}

\begin{dropdown}{_Proof_}
\begin{proof}
This is a proof. Prove things here. This should probably follow a theorem, lemma, or corollary.
This relates to \definitionref{named-label}.
\end{proof}
\end{dropdown}

\nbsp

\begin{algorithm}{Algo Name}
This is actually just like the lemma and theorem environments.

It is not tied to any of the standard latex pseudocode libraries.

This is intentional.
\end{algorithm}


If you need to use any other notation or symbols, you can do so in your submitted `.tex` file, but please avoid this where possible.

\theoremscripts

