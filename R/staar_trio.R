#' STAAR.Trio procedure using omnibus test on basis of a family-based design
#' @param genotype A 3n*p matrix for the original genotype data, in which n is the number of trios and p is the number of variants.
#' Each trio must consist of father, mother, and offspring (in this order). The genotypes must be coded as 0, 1, or 2. Missing genotypes are not allowed.
#' @param maf minor allele frequency
#' @param pos A numeric vector of length p for the position of p variants.
#' @param annotation_phred a data frame or matrix of functional annotation data
#' of dimension p*q (or a vector of a single annotation score with length p).
#' @param y_res A numeric vector of length n for the residual Y-Y_hat. Y_hat is the predicted value from the regression model in which the quantitative trait Y is regressed on the covariates.
#' If Y is dichotomous, you may treat Y as quantitative when applying the regression model.
#' @param xchr A logical value indicating whether the analysis is for the X chromosome. When xchr is TRUE, the analysis is for the X chromosome and sex is required.
#' When xchr is FALSE, the analysis is for the autosomes.
#' @param sex A numeric vector of length n for the sex of offspring. 0s indicate females and 1s indicate males. Sex is required when xchr is TRUE.
#' @param annotation_phred a data frame or matrix of functional annotation data
#' of dimension p*q (or a vector of a single annotation score with length p).
#' @param rare_maf_cutoff the cutoff of maximum minor allele frequency in
#' defining rare variants (default = 0.01).
#' @param rv_num_cutoff the cutoff of minimum number of variants of analyzing
#' a given variant-set (default = 2).
#' @return A list for the analysis results.
#' @importFrom stats pcauchy pnorm dbeta
#' @export

staar_trio <- function(genotype,maf,pos,annotation_phred=NULL,
                       rare_maf_cutoff=0.05,rv_num_cutoff=2,y_res=NA,xchr=FALSE,sex=NA
){

  if(class(genotype) != "matrix" && !(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix")){
    stop("genotype is not a matrix!")
  }

  if(dim(genotype)[2] == 1){
    stop(paste0("Number of rare variant in the set is less than 2!"))
  }



  n<-nrow(genotype)/3
  index_dad<-seq(1,(3*n),3) #index for parent1
  index_mom<-seq(2,(3*n),3) #index for parent2
  index_off<-seq(3,(3*n),3) #index for offspring



  annotation_phred <- as.data.frame(annotation_phred)
  if(dim(annotation_phred)[1] != 0 & dim(genotype)[2] != dim(annotation_phred)[1]){
    stop(paste0("Dimensions don't match for genotype and annotation!"))
  }

  MAF <- maf
  RV_label <- as.vector((MAF<rare_maf_cutoff)&(MAF>0))
  Geno_rare <- genotype[,RV_label]

  rm(genotype)
  gc()
  annotation_phred <- annotation_phred[RV_label,,drop=FALSE]


  if (!xchr){
    Z<-Geno_rare[index_off,,drop=F]-(Geno_rare[index_dad,,drop=F]+Geno_rare[index_mom,,drop=F])/2
  } else Z<-xcontribution(Geno_rare,sex)


  Z<-sweep(Z,1,y_res,'*')


  out<-list()

  out$Z<-Z


  if(sum(RV_label) >= rv_num_cutoff){
    #G <- as(Geno_rare,"dgCMatrix")
    MAF <- MAF[RV_label]
    rm(Geno_rare)
    #gc()
    pos = pos[RV_label]
    annotation_rank <- 1 - 10^(-annotation_phred/10)

    weight<-1/(sqrt(n*MAF*(1-MAF))) #MAF-based weights
    weight[which(weight==Inf)]<-0

    Z<-sweep(out$Z, 2, weight, '*')
    Zsq<-t(Z)%*%Z
    Zsum<-apply(Z,2,sum)
    ind<- 1:length(MAF)
    z<-fbat_z(Zsum,Zsq,ind) #FBAT burden z score for the ith window
    p_knockoff<-2*pnorm(-abs(z))

    ## beta(1,25)
    w_1 <- dbeta(MAF,1,25)

    Z<-sweep(out$Z, 2, w_1, '*')
    Zsq<-t(Z)%*%Z
    Zsum<-apply(Z,2,sum)
    z<-fbat_z(Zsum,Zsq,ind) #FBAT burden z score for the ith window
    p_w1<-2*pnorm(-abs(z))
    ## beta(1,1)
    w_2 <- dbeta(MAF,1,1)

    Z<-sweep(out$Z, 2, w_2, '*')
    Zsq<-t(Z)%*%Z
    Zsum<-apply(Z,2,sum)
    z<-fbat_z(Zsum,Zsq,ind) #FBAT burden z score for the ith window
    p_w2<-2*pnorm(-abs(z))


    ## Burden
    w_b <- sweep(annotation_rank,1,weight,"*")
    w_B_1 <- sweep(annotation_rank,1,w_1,"*")

    w_B_2 <- sweep(annotation_rank,1,w_2,"*")

    p.fbat <- c()
    p1.fbat<-c()

    p2.fbat<-c()
    # AW FBAT burden
    for (k in seq(dim(annotation_rank)[2])){
      Z<-sweep(out$Z, 2, as.vector(w_b[,k]), '*')
      Zsq<-t(Z)%*%Z
      Zsum<-apply(Z,2,sum)
      ind<- 1:length(MAF)
      z<-fbat_z(Zsum,Zsq,ind) #FBAT burden z score for the ith window
      p.fbat[k]<-2*pnorm(-abs(z))

      # w1 & w2
      Z1<-sweep(out$Z, 2, as.vector(w_B_1[,k]), '*')
      Zsq1<-t(Z1)%*%Z1
      Zsum1<-apply(Z1,2,sum) #
      Z2<-sweep(out$Z, 2, as.vector(w_B_2[,k]), '*')
      Zsq2<-t(Z2)%*%Z2
      Zsum2<-apply(Z2,2,sum) #
      z1<-fbat_z(Zsum1,Zsq1,ind) #FBAT burden z score
      p1.fbat[k]<-2*pnorm(-abs(z1)) #FBAT burden p-value
      z2<-fbat_z(Zsum2,Zsq2,ind)
      p2.fbat[k]<-2*pnorm(-abs(z2))

    }

    p_result <- list()
    p_result$p.fbat <- p.fbat
    p_result$p1.fbat <- p1.fbat
    p_result$p2.fbat <- p2.fbat
    p_result$p_0.5_0.5 <- p_knockoff
    p_result$p_w1 <- p_w1
    p_result$p_w2 <- p_w2

    return(p_result)}else{
      stop(paste0("Number of rare variant in the set is less than ",rv_num_cutoff,"!"))
    }
}

xcontribution <- function(dat,sex=NA){
  #sex: gender info for the offspring, 0 for females and 1 for males
  #Only knockoffs of mothers are required
  n.row <- nrow(dat)
  dad <- dat[seq.int(1, n.row, 3),, drop=FALSE]
  mom <- dat[seq.int(2, n.row, 3),, drop=FALSE]
  kid <- dat[seq.int(3, n.row, 3),, drop=FALSE]

  Z<-array(0,dim=c(n.row/3,ncol(dat)))
  het <- (mom == 1)
  hethom <- het & (dad == 0)
  Z[hethom & (kid == 0)]<-(-0.5)
  Z[hethom & (kid == 1)]<-0.5
  hethom <- het & (dad == 1)
  ind111<-hethom & (kid == 1)
  Z[hethom & (kid == 0)]<-(-0.5)
  Z[sweep(ind111,1,sex==0,'&')]<-(-0.5)
  Z[sweep(ind111,1,sex==1,'&')]<-0.5
  Z[hethom & (kid == 2)]<-0.5

  return(Z)
}


fbat_z<-function(W,V,ind){
  if (length(ind)>0){
    s1<-sum(W[ind])
    s2<-sum(V[ind,ind])
    if (s2==0) s2<-1
    return (unname(s1/sqrt(s2)))
  } else return(NA)
}
