---
title: Response to the Reviewers
author: Jakob Richter
header-includes: |
  \usepackage{framed}
  \newlength{\leftbarwidth}
  \setlength{\leftbarwidth}{3pt}
  \newlength{\leftbarsep}
  \setlength{\leftbarsep}{-10pt}
  \newcommand*{\leftbarcolorcmd}{\color{leftbarcolor}}
  \colorlet{leftbarcolor}{gray}
  \renewenvironment{leftbar}{%
      \def\FrameCommand{{\leftbarcolorcmd{\vrule width \leftbarwidth\relax\hspace {\leftbarsep}}}}%
      \MakeFramed {\advance \hsize -\width \FrameRestore }%
    \noindent
  }{%
      \endMakeFramed
  }
  \renewenvironment{quote}
  {\vspace{5pt} \leftbar\list{}{\rightmargin\leftmargin}%
  \item\vspace{-5pt}\relax}
    {\endlist\endleftbar}
---

# Second response to the reviews for "Improving Adaptive Seamless Designs through Bayesian optimization"

We would like to thank all reviewers again for their valuable and constructive comments on the paper.
We carried out a new simulation to address the comments by Reviewer 2 and included the results in the appendix.

## Reviewer: 1

### Comments to the Author

Thanks for addressing/responding to all my comments

## Reviewer: 2

### Comments to the Author
Thanks authors for the response to the comments. However, some of them are not addressed appropriately here. Please consider them carefully in the next round revision.
For page 7 table 1 effect size scenario, although it is common that early come is overoptimistic and higher, however, it is also often in practice that there are delayed treatment effect (e.g. in the immunotherapy clinical trial setting) where stage 1 effect are smaller than stage 2. I would strongly recommend authors to add difference simulation setting.

> We carried out new simulations to address this comment. 
> In the desired scenario, titled `paper_rev`, effect sizes for early and final are flipped compared to the standard scenario `paper`, so that the lower values occur in the early stage.
> The results are shown in the Appendix.

Similarly for the selection of “espilon” and “threshold” in page 8 table 2, there should be more practical evaluation where they are commonly based on the clinical meaningful cut, besides as optimization parameters. Suggest to add relevant discussion here to facilitate more practical usage of proposed methods.

> In this paper we consider an optimization problem maximizing the power given a total sample size.
> If, for example, the task was to select all treatments with clinically relevant effect sizes, then clearly the selection parameter needs to be informed by clinical reasoning.
> We added a respective statement to the discussion.

Page 8 table 3, authors should add more details or rationale on the how parameters range of search space is determined.

> We included a statement that “epsilon” and “threshold” are obtained through preliminary studies.
> Furthermore, we added a Figure S4 in the Appendix, which demonstrates that the ranges are chosen in a sensible way.

### Minor comments:
There is citation format typo in page 3: “as implemented in the asd package (?)”

> This occurred only in the diff-version, which will not be published.
