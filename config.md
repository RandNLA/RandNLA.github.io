<!--
Add here global page variables to use throughout your website.
-->
+++
author = "Septimia Zenobia"
mintoclevel = 2

# Add here files or directories that should be ignored by Franklin, otherwise
# these files might be copied and, if markdown, processed by Franklin which
# you might not want. Indicate directories by ending the name with a `/`.
# Base files such as LICENSE.md and README.md are ignored by default.
ignore = ["node_modules/"]

# RSS (the website_{title, descr, url} must be defined to get RSS)
generate_rss = true
website_title = "RandNLA Proof Wiki"
website_descr = "An incomplete list of the shortest, cleanest, and most important known proofs in the RandNLA literature"
website_url   = "https://randnla.github.io/"
+++

<!--
Add here global latex commands to use throughout your pages.
-->

\newcommand{\nbsp}{~~~&nbsp;~~~}
\newcommand{\highlight}[1]{\html{<span class="highlight">}!#1\html{</span>}}
\newcommand{\html}[1]{~~~!#1~~~}

\newenvironment{figure}[1]{
  \html{<figure>}
}{
  \html{<figcaption>#1</figcaption></figure>}
}

\newcommand{\includegraphics}[1]{![](!#1)}
