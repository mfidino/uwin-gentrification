# Notes on:

```
@article{mokany2011combining,
  title={Combining $\alpha$-and $\beta$-diversity models to fill gaps in our knowledge of biodiversity},
  author={Mokany, Karel and Harwood, Thomas D and Overton, Jacob McC and Barker, Gary M and Ferrier, Simon},
  journal={Ecology letters},
  volume={14},
  number={10},
  pages={1043--1051},
  year={2011},
  publisher={Wiley Online Library}
}
```

Their model is called `DynamicFOAM` (Dynamic framework for occurrence allocation in metacommunites).

1. Uses pair-wise compositional dissimilarity as the complement of sorensen's similarity coefficient (i.e., sorensens's dissimilarity `\beta_{ij} = 1 - (1c_{ij}/a_i+a_j`).

`a_i` is number of species as site `i`

`a_j` is number of species at site `j`

`c_{ij} is number of species in common between the two sites.`

**idea**: Even if we use `JAGS` we could calculate this with matrix math.

```R
zmat <- matrix(sample(0:1, 100, TRUE), ncol = 10, nrow = 10)

cij <- j %*% t(zmat)

# get richness
rich <- rowSums(zmat)

# matrix of same dimensions as cij
bij <- cij

for(i in 1:nrow(zmat)){
  for(k in 1:nrow(zmat)){
    bij[i,k] <- 1 - ((2 * cij[i,k]) / (rich[i] + rich[k])) 
  }
}
round(bij,2)

     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
 [1,] 0.25 0.33 0.25 0.20 0.33 0.33 0.27 0.45 0.23  0.09
 [2,] 0.78 0.60 0.56 0.64 0.60 0.40 0.67 0.67 0.57  0.50
 [3,] 0.50 0.78 0.75 0.40 0.78 0.56 0.45 0.64 0.54  0.45
 [4,] 0.80 0.64 1.00 0.67 0.64 0.27 0.38 0.38 0.47  0.54
 [5,] 0.78 0.40 0.56 0.64 0.40 0.60 0.50 0.50 0.43  0.67
 [6,] 0.56 0.60 0.33 0.27 0.20 0.20 0.50 0.33 0.29  0.33
 [7,] 0.27 0.67 0.45 0.38 0.17 0.33 0.29 0.43 0.25  0.29
 [8,] 0.27 0.50 0.45 0.38 0.50 0.67 0.43 0.29 0.25  0.14
 [9,] 0.54 0.57 0.54 0.33 0.57 0.43 0.38 0.25 0.22  0.25
[10,] 0.64 0.83 0.82 0.69 0.83 0.50 0.57 0.71 0.62  0.43
```

## Other requirements

1. \gamma diversity must be known (or estimated across all communites in a region). We can do with really easily with a multi-region model.


## Thoughts

- It is a permutation based approach that takes the input from different from different seperate alpha and beta diversity analyses. It's cool, but I'm not sure it's what we are looking for yet.

# Notes on:

```
@article{ferrier2007using,
  title={Using generalized dissimilarity modelling to analyse and predict patterns of beta diversity in regional biodiversity assessment},
  author={Ferrier, Simon and Manion, Glenn and Elith, Jane and Richardson, Karen},
  journal={Diversity and distributions},
  volume={13},
  number={3},
  pages={252--264},
  year={2007},
  publisher={Wiley Online Library}
}
```

This paper points out two approaaches that people use to model spatial patterns in biodiveristy:

1. 'predict, first assemble later': The modeled distributions of species are numerically classified to derive community types.

2. 'assemble first, predict later': communities are derived through the classificaiton of raw survey data, and then these are modeled.

I think we would take the first approach.

However, they suggest modeling emergent properties of biodiversity (e.g,. species richness).

Two ways to do that for beta-diversity:

1. the raw-data approach, partitions environmental / geographical components of beta diversity via canonical analysis (e.g., CCA). 

2. the 'distance' approach, dissimilarities in composition betwen pairs of sites are relative to environmental or geographical distances via matrix correlation or regression.

Second approach is better for analyzing variation among communities (which is what we want).

# Linear matrix regression

evaluates the correlation between two distance matrices.

They used Bray-Curtis dissimilarity:

`dij = 1 - (2A) / (2A + B + C)`

`A` = species common to both sites
`B` = number of species at site i
`C` = number of species at site j

Matrix regression can then be formulated as:

`d[i,j] = a0 + a1 *|x[i] - x[j|`

The linear model isn't ideal though, because the response variable is bounded between 0 and 1.

# The GLM version of matrix regression

They prefer the cloglog link by the looks of it. Inverse link function:

`u = 1 - exp(-eta)`

with binomial variance function

`V(U) u(1-u)/(si + sj)`

so response variable, it's complement, divided by the sum of species at a site.

They suggest using I-splines to estimate non-stationarity in rates of compositional turnover, and constrain the regression coefficient to be positive.

This model looks really cool, I should explore how to incorporate it into a Bayesian model.

## Notes on:

```
@article{woolley2017characterising,
  title={Characterising uncertainty in generalised dissimilarity models},
  author={Woolley, Skipton NC and Foster, Scott D and O'Hara, Timothy D and Wintle, Brendan A and Dunstan, Piers K},
  journal={Methods in Ecology and Evolution},
  volume={8},
  number={8},
  pages={985--995},
  year={2017},
  publisher={Wiley Online Library}
}
```

# Their model

y_ij ~ Binomial(n_ij, pi_ij)

y_ij = number of species not shared between sites i and j

n_ij = Total number of unique species at sites i and j.

pi_ij = expected dissimilarity

uses I-splines and constrains parameters to be positive, because larger distance should mean more dissimilar.
