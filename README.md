
# OrdCD

OrdCD is an R package (available on CRAN) for discovering causality from observational ordinal categorical data, developed and maintained by Yang Ni at Texas A&M University.

The package can also be downloaded at https://web.stat.tamu.edu/~yni/files/OrdCD_1.1.2.tar.gz.


#### Reference:  

Ni, Y., & Mallick, B. (2022). [Ordinal Causal Discovery.](https://proceedings.mlr.press/v180/ni22a/ni22a.pdf) *Proceedings of the 38th Conference on Uncertainty in Artificial Intelligence*, (UAI 2022), PMLR 180:1530–1540.


## A small simulation example
This example generates a network with 5 nodes, 3 categoreis, and 1000 observations. The true graph is a Markov chain.
```{r dataload,echo=TRUE,  warning=FALSE, message=TRUE }
set.seed(2020)
n=1000 #sample size
q=5 #number of nodes
y = u = matrix(0,n,q)
u[,1] = 4*rnorm(n)
y[,1] = (u[,1]>1) + (u[,1]>2)
for (j in 2:q){
  u[,j] = 2*y[,j-1] + rnorm(n)
  y[,j]=(u[,j]>1) + (u[,j]>2)
}
A=matrix(0,q,q) #true DAG adjacency matrix
A[2,1]=A[3,2]=A[4,3]=A[5,4]=1
y=as.data.frame(y)
for (j in 1:q){
  y[,j]=as.factor(y[,j])
}

G=OrdCD(y) #estimated DAG adjacency matrix
print(A)
print(G)
```

