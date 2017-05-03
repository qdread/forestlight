# Power law cutoffs

QDR, 2 May 2017

Clearly the BCI 50 hectare plots do not follow a power law at large sizes. The thinking is that there is a cutoff. Trees below the cutoff diameter are subject to density-dependent mortality, but trees above the cutoff diameter are subject to random external mortality. I believe the prediction should be that gap trees will have a higher cutoff than shade trees, since they are subject to density-dependent mortality up to a very big size. The eyeball test shows that the cutoff should be around 10<sup>1.5</sup>, or about 30 cm. 

The Caroline Farrior paper uses a piecewise function to find the cutoff. I have not yet figured out how to get that to work, but I did implement a very similar method. I fit a "power law with exponential cutoff" distribution, which is described in Clauset et al. I got the likelihood function from Cosma Chalizi's website, `http://www.stat.cmu.edu/~cshalizi/`. This function has three parameters: the scaling parameter `$alpha$`, the lower bound `$x_{min}$`, and the cutoff `$L$`. The function, which is a Pareto multiplied by an exponential, behaves like a Pareto power law below the cutoff, and like an exponential above.  

First, I figured out that the function is not available in the `poweRlaw` package. So I had to write the likelihood estimator out by hand. I wrote the basic Pareto one to confirm that it gave the same scaling exponent as the function from `poweRlaw`, which it did. I wrote the likelihood function for the Pareto with cutoff as well. The result, as we might have predicted, is that the estimated cutoff parameter increases going from shade-tolerant to intermediate to gap.

The columns of the table below are as follows:

- **slope1** is the density scaling parameter without a cutoff, for all trees. This is the same one we got earlier.
- **slope2** is the density scaling parameter for the portion of the trees below the cutoff in the Pareto-with-cutoff model. It is shallower because it isn't being dragged down by the big trees' random mortality anymore.
- **cutoff** is the diameter at which the exponential takes over from the power law (measured in units of cm). The cutoff is larger for the gap trees.
- **logcutoff** is just the `$log_{10}$` of the cutoff. It roughly corresponds to the eyeballed value.

| guild          | slope1 | slope2 | cutoff | logcutoff |
|----------------|--------|--------|--------|-----------|
| all            | 1.884  | 1.462  | 23.129 | 1.364     |
| shade-tolerant | 1.956  | 1.300  | 10.892 | 1.037     |
| intermediate   | 1.758  | 1.205  | 20.866 | 1.319     |
| gap            | 1.995  | 1.802  | 55.042 | 1.741     |

## Slope confidence intervals

Added 3 May 2017.

I did bootstrap sampling 99 times on the original datasets and then ran the MLE on each of those datasets to get a confidence interval around the slope estimates. I think it might not be a conservative enough method since the confidence intervals are very small&mdash;they don't even show up on the figure below except for the one around the breakpoint for gap trees. But I think the cutoffs are so different that any method would show that they are significantly different.

![slope figure](file:///C:\\Users\\Q\\Dropbox\\projects\\forestlight\\bootstrapci_cutoffs.png)
