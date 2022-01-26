# Principal Component Analysis

Principal Component Analysis (PCA) is one of the most popular techniques to perform "dimensionality reduction" of complex data sets. If we see the data with many variables as points in a high-dimensional space, we can compute new variables as linear combinations of the original ones and represent each data point as a set of coordinates in the new variables. In this way, we can project large-dimensional data sets onto low-dimensional spaces and lose the least information about the data.

## Motivation: simplifying complex data

Suppose we have a data set with $n$ variables and $m$ observations of each (typically, with $n \gg m$), in which the $m$ rows are observations and the $n$ columns are the variables. Each row of this matrix defines a point in the Euclidean space $\mathbb R^n$. Many biological data sets, e.g. gene expression {numref}`fig-micro-array`, RNAseq, medical imaging, can contain thousands or more variables, which poses major challenges both for visualization and computational tasks. PCA provides the best representation of such a data set in terms of a smaller set of variables, while capturing as much variance as possible. 

```{figure} figs/micro_array.jpg
---
name: fig-micro-array
---
Image of a microarray plate, [source](http://exploreable.files.wordpress.com/2011/04/array.jpg). Here each dot is a different variable (different gene) and this image in just one set of observations that will be placed into a row of the data matrix.
```

The intuition behind finding these new collective variables rests on the fact that the original variables have relationships. This is typically measured using covariance, which quantified how much a pair variables tends to move in the same direction (positive covariance) or in opposite directions (negative covariance). If two variables are tighlty coupled, one can replace the two measurements with one, which will describe how much the two of them are deviating in some collective way. 

It is helpful to think of this geometrically: if the variables are related, the scatterplot of observed data points will have a shape. The goal of PCA is to find directions in the $n$-dimensional space of observations that best match the shape of the data cloud.


## PCA algorithm

We start with a data set $X$ in the form of a $m$ by $n$ matrix. The first step is to decide which are the variables and which are the observations. For example, in the case of the microarray experiment, it usually makes sense to consider different genes the variables, and to use principal components to see which genes tend to be expressed together with others.


The second step is to compute the variance-covariance matrix of the $N$ variables. 

```{admonition} Definition
The *variance-covariance* matrix $C$ of a data set $X$ with $n$ variables $x_i$ and $m$ observations is an $n$ by $n$ matrix that contains pairwise variances between all $n$ variables, so that its $i$, $j$ element is:

$$C_{i,j} = Cov(X_i,X_j)$$

```
The third step is to diagonalize (find the eigenvalues and eigenvectors) of the covariance matrix $C$. The eigenvectors are the principal components of the $n$ variables in the data set, representing linear combinations of the variables that best fit the data. Diagonalizing an $n$ by $n$ matrix results in $n$ eigenvectors, so in order to simplify the description one needs to choose the most significant ones. This is accomplished by choosing a subset of $k$ principal components with the largest eigenvalues. Here are the steps of principal component analysis (PCA):

```{admonition} PCA algorithm
:class: tip
1. Obtain a dataset as a $m$ by $n$ matrix, with $n$ variables and $m$ observations
2. Compute covariances for variable i and variable j, put them in the covariance matrix $C$
3. Compute the eigenvalues and eigenvectors (principal components) of the matrix $C$
4. Order the principal component by size of eigenvalues from largest to smallest and select a few as the new coordinates
```

## Optization by explained variance

The reason that we order the PCs by their eigenvalues is that they measure the amount of variance captured by each principal component. In that, they are equivalent to the coefficient of determination $r^2$  in linear regression. The sum of all the eigenvalues is equal to the total variance of all the variables:

$$ \sum_i \lambda_i = \sum Var(X_i)$$

and the fraction of variance captured by the a principal component is:

$$ Var(PC_i) = \frac{\lambda_i}{\sum_i \lambda_i}$$

The theory behind this rests on some relatively sophisticated linear algebra, in particular what is called the singular value decomposition (SVD) and the Eckart-Young Mirsky theorem. Here is a nice video by Gilbert Strang that explains this:[Strang lecture](https://www.youtube.com/watch?v=Y4f7K9XF04k)


## Dimensionality reduction

After sorting the principal components and selecting $k$ largest eigenvalues, we are ready to simplify the data. This means that we can express a data set of $n$ variables in terms of the coordinate set of $k$ principal components. In order to express the data set in this new system of coordinates, we compute the projection coefficients for each measurement onto a give principal component. Suppose that $Y$ is a set of measurements of $N$ variables (e.g. genes) and $P_i$ is the $i$-the principal component. Then the projection coefficient of $Y$ onto $P_i$ is the dot product of the two vectors (both of length $N$) divided by the squared norm (length) of the PC:

$$ c_i = \frac{\langle Y, P_i\rangle}{|| P_i ||^2} $$

If the eigenvectors are normalized prior to the computation (as they are by most computational packages), then the projection coefficient is just the dot product. Then the coefficients can be obtained for all of the measurements in the data set $X$ ($m$ by $n$) by multiplying it by the matrix $P$ containing the first $k$ eigenvectors (principal components), which has $n$ rows and $k$ columns. The result is an $m$ by $k$ matrix $C$ containing $k$ coefficients for each of the $m$ measurements:

$$
C = D \times P
$$


Here is the outline for the transformation:
```{admonition} Dimensionality reduction
:class: tip
1. Subtract the mean of each observation from the data matrix (if it has $M$ observations in rows and $N$ variables as columns, subtract the mean of each row from it)
2. Compute the projection coefficients $C = D \times P$ for each measurement and each of the $k$ principal components
3. Plot or otherwise display these coefficients as coordinates in the new vector system of the $k$ PCs. This can be used to cluster or otherwise find patterns in the observations.
```


The entire data set can be expressed in a low-dimensional setting, for instance plotted in the plane of two principal components with coordinates $(c_{i,1}, c_{i,2})$ for each data measurement $i$. This is often useful for clustering, or grouping experimental conditions based on the similarity of their principal component representations. Biologists frequently do this with complex data sets, for example grouping different cell lines together by their gene expression profiles.


