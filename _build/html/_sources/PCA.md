# Principal Component Analysis
Oftentimes in modern biology you will encounter a very large data set, with many variables and even more observations. In order to simplify it and extract some interesting information, one can find its principal components. 

## Motivation: simplifying data


## PCA algorithm

To find the principal components, one needs to start with a data set in the form of a $M$ by $N$ matrix. For example, this could be $M$ gene expression profiles measured under $N$ different conditions in a microarray experiment, as in {numref}`fig-micro-array`.The first step is to decide which are the variables and which are the observations. In the case of the microarray experiment, it usually makes sense to consider different genes the variables, and to use principal components to see which genes tend to be expressed together with others.


```{figure} figs/micro_array.jpg
---
name: fig-micro-array
---
Image of a microarray plate, [source](http://exploreable.files.wordpress.com/2011/04/array.jpg)
```


The second step is to compute the variance-covariance matrix of the $N$ variables. 

```{admonition} Definition
The *variance-covariance* matrix $C$ of a data set with $N$ variables $X_i$ and $M$ observations is an $N$ by $N$ matrix that contains pairwise variances between all $N$ variables, so that its $i$, $j$ element is:
$$ C_{i,j} = Cov(X_i,X_j)$$
```


The third step is to diagonalize (find the eigenvalues and eigenvectors) of the covariance matrix $C$. The eigenvectors are the principal components of the $N$ variables in the data set, representing linear combinations of the variables that best fit the data. Diagonalizing an $N$ by $N$ matrix results in $N$ eigenvectors, so in order to simplify the description one needs to choose the most significant ones. This is accomplished by choosing a subset of $k$ (usually 2-5) principal components with the largest eigenvalues. Here are the steps of principal component analysis (PCA):

```{admonition} PCA algorithm
:class: tip
1. Obtain dataset as a $M$ by $N$ matrix, with $N$ variables and $M$ observations
2. Compute covariances for variable i and variable j, put them in the covariance matrix $Cov$
3. Compute the eigenvalues and eigenvectors (principal components) of the matrix $Cov$
4. Order the principal component by size of eigenvalues from largest to smallest and select a few as the new coordinates
```

## Dimensionality reduction

The reason that we order the PCs by their eigenvalues is that they measure the amount of variance captures by each principal component. In that, they are equivalent to the coefficient of determination $r^2$  in linear regression. The sum of all the eigenvalues is equal to the total variance of all the variables:

$$ \sum_i \lambda_i = \sum Var(X_i)$$

and the fraction of variance captured by the a principal component is:

$$ Var(PC_i) = \frac{\lambda_i}{ \sum_i \lambda_i}$$

After sorting the principal components and selecting 2-5 with the largest eigenvalues, what do we do? First, as stated, they are useful for reducing the dimensionality of a data set. This means that one can choose to express a data set of $N$ variables in terms of the coordinate set of $k$ principal components. In order to express the data set in this new system of coordinates, we compute the projection coefficients for each measurement onto a give principal component. Suppose that $Y$ is a set of measurements of $N$ variables (e.g. genes) and $P_i$ is the $i$-the principal component. Then the projection coefficient of $Y$ onto $P_i$ is the dot product of the two vectors (both of length $N$) divided by the squared norm (length) of the PC:

$$ c_i = \frac{\langle Y, P_i\rangle}{|| P_i ||^2} $$

If  the eigenvectors are  normalized prior to the computation (as they are by most computational packages), then the projection coefficient is just the dot product. Then the coefficients can be obtained for all of the measurements in the data set $D$ ($M$ by $N$) by multiplying it by the matrix $P$ containing the first $k$ eigenvectors (principal components), which has $N$ rows and $k$ columns. The result is an $M$ by $k$ matrix $C$ containing $k$ coefficients for each of the $M$ measurements:

$$
C = D \times P
$$

Thus,  the entire data set can be expressed in a low-dimensional setting, for instance plotted in the plane of two principal components with coordinates $(c_{i,1} c_{i,2})$ for each data measurement $i$. This is often useful for clustering, or grouping experimental conditions based on the similarity of their principal component representations. Biologists frequently do this with complex data sets, for example grouping different cell lines together by their gene expression profiles.


Here is the outline for the transformation:
```{admonition} Dimensionality reduction
:class: tip
1. Subtract the mean of each observation from the data matrix (if it has $M$ observations in rows and $N$ variables as columns, subtract the mean of each row from it)
2. Compute the projection coefficients $C = D \times P $ for each measurement and each of the $k$ principal components
3. Plot or otherwise display these coefficients as coordinates in the new vector system of the $k$ PCs. This can be used to cluster or otherwise find patterns in the observations.
```



