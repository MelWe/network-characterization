# -------------------------------------------------- #
# Sample script for computing Forman-Ricci curvature
#
# Author: Melanie Weber (mw25@math.princeton.edu)
#
# Please cite:
#
# @article{weber2017characterizing,
#    title={Characterizing Complex Networks with Forman-Ricci Curvature and Associated Geometric Flows},
#    author={Weber, Melanie and Saucan, Emil and Jost, J{\"u}rgen},
#        journal={Journal of Complex Networks},
#        volume={5},
#        number={4},
#        pages={527--550},
#        year={2017},
#        publisher={Oxford University Press}
# }
#
# @article{weber2016forman,
# title={Forman-Ricci flow for change detection in large dynamics data sets},
# author={Weber, Melanie and Jost, J{\"u}rgen and Saucan, Emil},
#    journal={Axioms},
#    volume={5},
#    number={4},
#    pages={doi--10},
#    year={2016}
# }
#
# -------------------------------------------------- #

library(compiler)
library(parallel)

# --- functions --- #

FR_curve <- function(e,edges,list,w_e,w_n){
    ## Ricci-curvature (weighted)
    # e: edge index
    # edges: edgelist
    # w_e: edge weights
    # w_n: node weights
    
    v1 <- as.numeric(edges[e,1])
    v2 <- as.numeric(edges[e,2])
    
    i1 <- which(list==v1)
    i2 <- which(list==v2)
    
    e1 <- which(as.numeric(edges[,1]) == v1) # in v1
    e2 <- which(as.numeric(edges[,2]) == v2) # out v2
    
    sum1 <- as.numeric(0)
    sum2 <- as.numeric(0)
    
    if(length(e1)>0){
        for(i in 1:length(e1)){
            sum1 <- as.numeric(sum1) + as.numeric(w_n[i1]/(sqrt(abs(w_e[e]*w_e[as.numeric(e1[i])]))))
        }
    }
    
    if(length(e2)>0){
        for(i in 1:length(e2)){
            sum2 <- as.numeric(sum2) + as.numeric(w_n[i2]/(sqrt(abs(w_e[e]*w_e[as.numeric(e2[i])]))))
        }
    }
    
    ric <- (w_n[i1] + w_n[i2]) - w_e[e]*(as.numeric(sum1)+as.numeric(sum2))
    
    return(ric)
}

FR <- cmpfun(FR_curve)

# --- main --- #

# input/ adjust net
data <- read.csv(file="data.csv", header=TRUE, sep=",") # read in edge list
list <- append(data[,1],data[,2])
list <- list[which(!duplicated(list))]
weight_n <- rep(1,length(list)) # asign node weights; here: unweighted | alternative: weighted node-degrees
weight_e <- data[,3] # asign edge weights

# FR-curvature weighted
# given: {data,list,weight_e,weight_n}
ind <- seq(1,length(weight_e),1)
FR_w <- mclapply(ind,FR,data,list,weight_e,weight_n,mc.cores=15L) # parallel computation with 15 cores
FR_w <- as.vector(unlist(FR_w))
FR_w2 <- abs(FR_w)/max(abs(FR_w))

# curvature plot (histogram)
png(file="FR.png",width=640,height=480)
hist(FR_w2,xlab="FR curvature",ylab="distribution",
breaks=50,xlim=c(0,1),
main=" ",
col="gray64",prob=T,
cex.lab=1.5,cex.axis=1.5,bty='n')
dev.off()


