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

# Response to the reviews for "Improving Adaptive Seamless Designs through Bayesian optimization"

We would like to thank all reviewers for their valuable and constructive comments on the paper.
We especially thank Reviewer 1 for the thoughtful feedback that helped us to clarify the definition of the total number of samples and the power.

## Reviewer: 1

### Comments to the Author

This is a well organised and written manuscript, describing a very thorough investigation of a reasonably interesting problem. Computational aspects of complex adaptive clinical trial designs have all too often been neglected. The speed gain proposed here is substantial, with no apparent loss in quality. The findings could help improve the uptake of adaptive seamless designs. I only have a few relatively minor suggestions that I hope the authors will find helpful.

### General

1) Power is used as a key performance metric, yet it lacks an explicit definition. The problem, of course, is that power is slightly more complicated in a multi-stage design with multiple arms (and therefore multiple hypotheses to be tested) than it is in a simple comparison of two groups. Given the scenario described in the motivating example, I would have assumed that power was defined as the probability of at least one (any) treatment group being demonstrated to be ‘better’ than control i.e. a ‘disjunctive’ or ‘minimal’ power definition (cf. Senn & Bretz 2007). However, on p. 10 it is stated that “the vector of treatment numbers for determining power counts rejections of 3 or 4 hypotheses”, implying a power definition more similar to ‘conjunctive’ or ‘maximal’ power.

Reference: Senn S, Bretz F (2007) Power and sample size when multiple endpoints are considered. Pharmaceutical Statistics, 6(3), 161-170.

> As in Friede et. al. (2020), where our motivating example was introduced, we count the rejections of any of the hypotheses 3 or 4 for calculating the power.
> We now clearly state this in the paper and also clarify that there are different options, and it is up to the user to determine the method of how the power is calculated.

2) There is some inconsistent notation:

a) It seems that n_treat is varyingly defined as the “total number of treatments” (p. 8, l. 29 and p. 9, l. 20), “total number of patients” (p. 8, l. 48 and p. 9, l. 14) and “allowed treatment number” (p. 12, l. 50 and p. 12, l. 51).

> We improved our wording, so that the number of treatment arms cannot longer be confused with the individual number of treatments of each patient.
> Therefore, the total sample size (formerly “total number of treatments”) is denoted as n_total (formerly n_treat).
> Since the patients are different in both stages (i.e. each patient is treated only once) the total number of patients is the same as the total sample size.
> However, we only use the phrase "total sample size" to avoid confusions.

b) On p. 3 k is introduced as the running index for treatments (k = 1, …, K) but later on p. 9 k is the resolution of the search grid.

> Thanks for catching that mistake. It is now fixed, we use l for the resolution of the search grid.

### Specific

p. 2: For the benefit of readers less familiar with the term, the authors should explain at least briefly what they mean by “expensive black-box” optimisation problems.

> We added a short explanation.

p. 3 ff.: This may be implicitly clear but it should also be said explicitly that in the motivating example and throughout the remainder of the paper a normality assumption is made for all outcome data.

> There is no normality assumption for the individual outcome data but for the test statistics Z_i,j a multivariate normal distribution is assumed which will hold approximately if the sample size is sufficiently large. 
> We clarified that the test statistics are multivariate normal distributed at least approximately at the very end of chapter 3.

p. 9: How were the values for the resolution of Grid and Grid Small (k=25 and k=7) selected?

> We added an explanation and referenced previous studies that used similar values.

p. 10: Similarly, how was the number of randomly sampled points (16) for the initial design for BO chosen?

> We added an explanation that often 10% of the budget is used for the initial design, together with a reference.

p. 10: In section 5.4 ten additional stochastic simulations are mentions, which seems in conflict with the 20 additional stochastic simulations mentioned in other parts of the manuscript.

> Thanks, for spotting that mistake. Indeed, we use 20 stochastic repetitions.

p. 10: It should be clarified whether the test level of 0.025 is one-sided.

> We use a one-sided test and added this information in the table.

p. 11: Whilst I appreciate that the focus is on relative improvements in runtime, the authors should at least briefly describe the hardware setup with which these runtimes were achieved, to provide some context for the absolute runtimes.

> We added the requested information.

p. 12: It is stated that: “Only for the scenarios (effect: paper2, n_treat = 1000) and (effect: sigmoid, n_treat = 1000) Grid Small has a small advantage over BO and Grid.” I think it should be added that Grid Small also has an advantage over BO for (effect: paper, n_treat = 1000) and certainly for (effect: sigmoid, n_treat = 2000).

> Thanks for bringing that up. We revised the list of cases where Grid Small shows a notable advantage.
> An advantage is not notable if the median of both methods lies within the box of the other method.
> Therefore, we do not see a notable advantage of Grid Small for the configuration (effect: paper, n_treat = 1000).

p. 12: It says: “For example, for scenario effect: sigmoid, n_treat = 1000 Grid is superior to BO, but the difference is smaller than 0.25%.” This isn’t true for n_treat = 1000 where the difference between Grid and BO is around 2%, but rather for n_treat = 2000.

> Thank you for spotting that mistake. This was changed.

p. 12: Apologies but I’m struggling to understand the following sentence: “As intuitively anticipated, in almost all scenarios the worst results were obtained if all available treatments are used in stage 1 and none in stage 2.” How can none of the treatments be used in stage 2 when on p. 8 it is specified that the number of treatments in stage 2 must always be between 2 and 5?

> If r=1, then all samples of the total number of allowed samples are evaluated in stage 1. 
> Therefore no more treatments can be evaluated in stage 2.
> We revised the wording to make this clearer.

p. 13: It is stated that: “Note, that the selection strategies thresh and all never selected powerful designs.” Why is that? Presumably, the former is true because the range of numerical values used for the threshold was chosen too narrow, whilst the latter shouldn’t be surprising at all, as it is the only non-adaptive design in the mix.

> In the analysis of the results we saw that the optimal threshold values are not close to the border.
> This indicates that the ranges are wide enough.
> We assume that the threshold selection is disadvantageous because it looks at the absolute values, which is not that adaptive even with optimization.

p. 14: What are the x-values referred to in the caption of Figure 4? No x is introduced in the notation, and Figure 4 itself plots y-values against r-values.

> Thanks for spotting this mistake. We changed the caption, replacing "x" by "configuration".

p. 16: The citation of the technical report by Bischl et al. should include a URL if possible.

> We added the arxiv id.

### Typos

> We corrected all spotted typos. Thank you.

## Reviewer: 2

### Comments to the Author

The authors propose to use Bayesian optimization (BO) to improve the efficiency of the design selection process in clinical trials. A set of parameters to be chose for the design optimization based on the power. The idea seem novel and could be a useful approach. However, there are several major issues and limitations of proposed methods: for example, power is often only part of measure of choosing design, other factors like study duration, and # of pts are also important. Also, the black-box function relay too much of the various parametric assumption of parameters, which need to be carefully decided to be more related to true clinical trial needs. Overall, I have the following comments for the authors.

> We agree that the power is only one measure for choosing a design.
> Therefore, the idea behind our method is to quickly asses which power one can achieve given the other fixed factors such as total sample size (and therefore the number of patients).
> Quickly obtaining the best possible performance for a specific design choice will also help decision makers to decide how to choose the design.
> Furthermore, we added study duration as a possible alternative performance measure. 
> This is now discussed as a limitation in our conclusion section.
> Assuming certain parameters under the null hypothesis is a standard approach.
> However, we see that a more Bayesian approach is certainly interesting, as we already stated in the last paragraph of the discussion.

### Major comments:

In page 7 table 1, the effect size scenarios presented seem always assume 2nd stage effect size are higher than stage 1? The authors should also evaluate the vice visa situation, and consistent effect size case in simulation.

> In the table the effects for the second stage (final) are indeed smaller than for the fist stage (early).
> This is a realistic scenario, as often the early outcome is overoptimistic and therefore higher.

In page 8 table 2, the selection of “espilon” and “threshold” used in the arm selection should be based on the clinical meaningful cut, I am not sure how practical that they are part of parameters to be optimized? The authors should also evaluate with different setting.

> Given a simulation setting we optimize the parameters epsilon (eps) and tau (thresh). This is independent of the application but of course the clinical application can give hints which simulation settings are adequate.

In section 5.2, the authors keep the total number of pts fixed, which seem not appropriate in practice, as overall power relying several factors, e.g. # of pts, the interim selection rule, stage 1 and stage 2 ratio. Often # of pts is also key goal of optimal design evaluation. I will strongly suggest # of pts should be a parameter in the black-box function, instead of fixing it to a constant value.

> If we included the number of patients in the optimization, the optimizer would propose to choose the maximal number of patients, because this gives the maximal performance.
> This is also indicated in our results, as we tried different values of n_total which equals the number of patients. 
> Here, a higher number always resulted in a higher power.
> As mentioned in the paper, the idea is to run this optimization for different sample sizes (denoted by n_total), and then to decide which sample size is required.
> This brings us to the future goal that we mentioned in the discussion: We can include n_total as another parameter and outcome in a multi-criteria optimization setting.

In page 10 table 4, the number of simulation iteration is only 1000, which does not seem sufficient for complex optimization methods proposed by authors. I would suggest at least 5000 should be needed.

> The simulation iterations are a specific parameter of the black-box function we optimize.
> The higher this value, the more reliable is the estimator of the power of the design.
> Our optimization is run for 100 iterations which is sufficient to converge as Figure 5 shows.
> Therefore, we do not see the need to run a larger number of iterations.
> We added a corresponding statement in section 5.7.

Page 11 table 5, the authors present the average run time in hour, not sure how these information is helpful, as this depend on the what computing systems were used in the simulation. Please either add more information.

> Thank you. We added the hardware setup that was used.

Page 13, BO has issue to find optimal value when one parameters closer to the borders of search space. This is concerning and confirm the limitation of using too many assumptions for optimization, e.g. appropriate parameters range. I am not sure taking log-transform could solve this problem, could authors elaborate more details on this.  Additional simulation is necessary to confirm this.   

> As we mentioned in the final discussion section of the paper, we see this as a topic for a future benchmark to further improve the performance of BO.
> The aim of this paper was to benchmark BO in a basic setting without the prior knowledge that small values of **r** are of higher importance.

Given authors used COPD trial as motivating example, I would strongly suggest to add section of case study for the illustration of proposed methods.

> In fact the effect set "paper" is based on a case study where the effect sets are estimated on a real dataset (Friede et al, 2020). 
> In this paper we do not propose a method for analyzing a real dataset, but a method to select the best design given a known theoretical situation.

### Minor comments:

Page 8 line 29,  n_{treat} should be total # of pts, instead of “total number of treatment”.

> Thanks, we fixed that mistake.
