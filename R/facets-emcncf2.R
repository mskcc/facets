#genotype mixture model using EM algorithm to call allele-specific copy number and cellular fraction
emcncf2=function(x,trace=FALSE,unif=FALSE,min.nhet=15,maxiter=10,difcf=0.05,maxk=5,eps=1e-3){  
  
  warning("emcncf2 occationally returns quirky copy number estimates due to the clonal cluster structure imposed. please use with caution.")
  
  jointseg=x$jointseg
  out=x$out
  dipLogR=x$dipLogR
  nX=x$nX
  seg=out
  
  jointseg=subset(jointseg,!is.na(jointseg$cnlr))  
  logR=jointseg$cnlr 
  #center logR at dipLogR
  logR.adj=logR-dipLogR
  
  logOR=jointseg$valor 
  logORvar=jointseg$lorvar 
  logOR2=logOR^2
  logOR2var=logOR2/logORvar
  logORvar.clust=by(logORvar,jointseg$segclust,function(x)mean(na.omit(x)))
  
  het=jointseg$het
  
  nmark=seg$num.mark
  segclust=seg$segclust
  nmark.clust=by(nmark,segclust,mean)
  
  cnlr.median.clust=by(seg$cnlr.median,segclust,function(x)mean(na.omit(x)))
  #mafR.clust=by(seg$mafR,segclust,function(x)mean(na.omit(x)))
  mafR.clust=by(seg$mafR.clust,segclust,function(x)mean(na.omit(x)))
  segs=rep(1:length(nmark),nmark)  
  nseg=length(nmark)
  nhet=seg$nhet
  nhet.clust=by(nhet,segclust,mean)
  chr=seg$chrom
  chr.clust=by(chr,segclust,mean)
  
  #if(nseg>500)stop("Likely hyper-segmented. Increase cval in procSample.")
  
  endseq=jointseg[cumsum(nmark),2]
  startseq=jointseg[c(1,cumsum(nmark)[-nrow(seg)]+1),2]
  seglen=(endseq-startseq)/1e6
  seglen.clust=by(seglen,segclust,mean)
  
  mafR=seg$mafR
  mafR[mafR<0]=0
  seglogr=seg$cnlr.median
  nclust=max(segclust)
  seglogr.clust=by(seglogr,segclust,mean)
  seglogr.clust.adj=as.vector(seglogr.clust)-dipLogR
  
  emflags=NULL
  
  var=var(jointseg$cnlr,na.rm=T)    
  if(var>0.6){
    logR=rep(seglogr,nmark)
    emflags=paste(emflags,"Noisy sample, Calls can be unreliable.",sep=" ")
  }
  
  #consider genotypes up to t=6, assume minor alelle is B, switch to simple moment estimates for high copy numbers (t>6) for computational efficiency
  genotype=c("0","A","AA","AB","AAB","AAA","AAAB","AABB","AAAA","AAAAB","AAABB","AAAAA","AAABBB","AAAABB","AAAAAB","AAAAAA")
  minor=c(0,0,0,1,1,0,1,2,0,1,2,0,3,2,1,0)
  major=c(0,1,2,1,2,3,3,2,4,4,3,5,3,4,5,6)
  t=ifelse(genotype=="0",0,nchar(genotype))
  
  ng=length(genotype)
  n=length(logR)  
  
  #diploid genome check 
  if(all(seg$cf[seg$chrom<nX]==1&seg$tcn[seg$chrom<nX]==2)|max(mafR.clust[seg$chrom<nX & seg$nhet>15], na.rm = T) < 0.05){
    rho=NA
    gamma=2
    out1=data.frame(seg[,1:9],start=startseq,end=endseq,cf.em=seg$cf,tcn.em=seg$tcn,lcn.em=seg$lcn)
    emflags=paste(emflags,"Insufficient information to estimate purity. Likely diplod or purity too low.",sep=" ")
    out=list(purity=rho,ploidy=gamma,dipLogR=dipLogR,cncf=out1, emflags=emflags)
    return(out)    
    stop("Insufficient information",call.=F)
  }
  
  
  #set intialize value of tumor purity from naive estimates
  rhov.lsd=seg$cf
  minor.lsd=seg$lcn
  t.lsd=seg$tcn
  major.lsd=t.lsd-minor.lsd
  nas=(is.na(major.lsd)|is.na(minor.lsd))
  
  homdel=which(major.lsd==0&major.lsd==0)
  genotype.lsd=rep(NA,nseg)
  
  a=lapply(1:sum(!nas),function(x)paste(rep("A",major.lsd[!nas][x]),collapse=""))
  b=lapply(1:sum(!nas),function(x)paste(rep("B",minor.lsd[!nas][x]),collapse=""))
  genotype.lsd[!nas]=unlist(lapply(1:sum(!nas),function(x)paste(a[[x]],b[[x]],sep="")))
  genotype.lsd[homdel]="0"
  which.geno.lsd=match(genotype.lsd,genotype)
  
  rhov.lsd[t.lsd==2&minor.lsd==1]=NA
  rhov.lsd[t.lsd==2&rhov.lsd==1]=NA
  rhov.lsd[chr>=nX&rhov.lsd==1]=NA
  
  naive=quantile(rhov.lsd,prob=0.75,na.rm=T)
  
  rhov.lsd.subset=rhov.lsd
  #avoid genotypes with identifiability issue
  rhov.lsd.subset[which.geno.lsd%in%c(3,5,7,10,11,14,15,NA)]=NA 
  rhov.lsd.subset[t.lsd>6]=NA
  rhov.lsd.subset[seglen<50]=NA
  loh=which(t.lsd>=1 & minor.lsd==0 & seglen>15) #use only LOH seg for initial estimate
  
   
  rho=NA
  if(length(loh)>2 &!all(is.na(rhov.lsd[loh]))){
    rho=max(by(rhov.lsd[loh],segclust[loh],function(x)mean(x,na.rm=T)),na.rm=T) #use LOH to estimate putity when available.
  }else{  
    if(length(na.omit(rhov.lsd.subset))>1)rho=max(by(rhov.lsd.subset,segclust,function(x)mean(x,na.rm=T)),na.rm=T)
  }
  
  if(is.na(rho)|rho<0.2)rho=naive
  
  rhov=rhov.lsd   
  rhov[is.na(rhov)]=rho
  #avoid initial value too low
  rhov[rhov<0.1]=rho
  #avoid 1
  rhov[rhov==1]=rho 
  rhov[nhet<min.nhet]=rho
  rhov[which.geno.lsd==4]=NA
  rhov=as.vector(by(rhov,segclust,function(x)mean(x,na.rm=T)))
  rhov0=rhov
  rhov0[rhov0>rho]=rho
  
  #intital single clone fit. set a single cf at purity for all segment clusters
  rhov=rep(rho,nclust)
  
  #initial value for genotype priror
  prior=matrix(1/ng,nrow=nclust,ncol=ng)
  posterior=rep(0.5,nclust)
  
  #initial value for sigma 
  sigma=rep(2,nclust)  
  
  #cold start for rho if set unif=TRUE
  if(unif){
    rhov=runif(nclust,0.3,0.8)
  }
  
  environment(onepass)=environment()
  
  #fit a single clonal cluster
  rho.clust=rep(1,nclust)
  rhov=rep(rho+0.05,nclust)
  cat("fitting 1 clonal cluster ...",'\n')
  em.out=onepass(x=x,trace=trace,unif=unif,maxiter=10,min.nhet=min.nhet,eps=eps,rho.clust=rho.clust,rho=rho,rhov=rhov,prior=prior, posterior=posterior, sigma=sigma)
  which.geno=em.out$which.geno
  posterior=em.out$posterior
  threshold=max(posterior[seglen.clust>10],na.rm=T)
  rhov=em.out$rhov
  rho=em.out$rho
  segclust.rhov=em.out$segclust.cf
  segclust.rhov[is.na(segclust.rhov)]=rho
  segclust.rhov[segclust.rhov>rho]=rho

  
  #Identify poor fits by posterior probability to refit as subclonal cluster. 
  cond1=(threshold-posterior)>0.05
  #The cf has to be sufficiently different to allow a subclonal cluster fit
  cond2=abs(segclust.rhov-rho)>difcf
  #Don't allow subclonal fit of tiny segments
  #cond3=seglen.clust>1
  cond3=(nmark.clust>50&nhet.clust>15)
  #Genotypes with identifiability issue don't allow subclonal option
  ub=c(10,15)
  cond4=(which.geno%in%ub)
  #Exclude sex chromosome
  cond5=chr.clust<nX
  #Exclude small changes
  cond6=(abs(seglogr.clust.adj)>0.15|mafR.clust>0.05)
  refit=which(cond2&cond3&!cond4&cond5&cond6)
  
  #recursively identify poor fits and fit additional subclonal clusters up to 4
  if(em.out$rho>0.4&any(refit)){
  nclone=2
  dif=refit
  difrho=abs(max(segclust.rhov[refit],na.rm=T)-rho)
  
  while(any(refit)&any(dif)&difrho>difcf&nclone<maxk){

  refit.old=refit
  rho.clust.start=rho.clust
  rho.clust.start[refit]=nclone
  
  #categories check to avoid overring
  #ncat=length(table(rho.clust))
  #if(max(rho.clust)>ncat)rho.clust[rho.clust>1]=rho.clust[rho.clust>1]-1
  
  rhos=by(rep(segclust.rhov,as.vector(nmark.clust)),rep(rho.clust.start,as.vector(nmark.clust)),function(x)mean(x,na.rm=T))
  rhov.start=rhos[rho.clust.start]
  rhov.start=as.vector(rhov.start)

  
  em.out=onepass(x=x,trace=trace,unif=unif,maxiter=maxiter,min.nhet=min.nhet,eps=eps,
                 rho.clust=rho.clust.start,rho=rho,rhov=rhov.start, prior=prior, posterior=posterior, sigma=sigma)
  
  posterior=em.out$posterior
  which.geno=em.out$which.geno
  threshold=quantile(posterior[seglen.clust>10],1,na.rm=T)
  
  rhov1=em.out$rhov
  rho1=em.out$rho
  rhos1=by(rep(rhov1,as.vector(nmark.clust)),rep(em.out$rho.clust,as.vector(nmark.clust)),function(x)mean(x,na.rm=T))
  
  segclust.rhov=em.out$segclust.cf
  segclust.rhov[is.na(segclust.rhov)]=rho1
  segclust.rhov[segclust.rhov>rho1]=rho1
  
  #set smaller threshold for posterior prob difference to have sufficient sensitivity
  cond1=(threshold-posterior)>0.05
  cond2=abs(segclust.rhov-rhov1)>difcf
  cond4=(which.geno%in%ub)
  refit=which(cond2&cond3&!cond4&cond5&cond6)
  dif=setdiff(refit.old,refit)
  difrho=min(abs(mean(rhov1[refit.old],na.rm=T)-rhos1[-nclone]))
  
  if(difrho>difcf){
    cat("fitting",nclone,"clonal clusters ...",'\n')
    rhov=rhov1
    rho=rho1
    rho.clust=em.out$rho.clust
    nclone=nclone+1
  }
  
  }
  }
  
  #rhov=em.out$rhov
  #rho=em.out$rho
  
  rhov.em=rhov[segclust]
  which.geno.em=which.geno[segclust]
  
  genotype.em=which.geno.em  
  genotype.em[which(!is.na(which.geno.em))]=paste("A",major[which.geno.em],"B",minor[which.geno.em],sep="")[which(!is.na(which.geno.em))]
  t.em=t[which.geno.em]
  major.em=major[which.geno.em]
  minor.em=minor[which.geno.em]
  
  clevel=unique(rhov.em)
  clevel=clevel[!is.na(clevel)&clevel!=1]
  corder=rank(-clevel)
  names(corder)=clevel
  
  clonal.cluster=corder[as.character(rhov.em)]
  clonal.cluster[is.na(rhov.em)]=NA
  clonal.cluster[which.geno.em==4]=1
  
  #clonal.cluster=rho.clust[segclust]
  
  #calculate ploidy
  gamma=(2^(-dipLogR)*(2*(1-rho)+2*rho)-2*(1-rho))/rho
  
  #hybrid: for high copy number (t>6), use moment estimates
  seglogr.adj=seg$cnlr.median-dipLogR
  idx=which(seglogr.adj>1.7*rho|is.na(which.geno.em))
  if(any(idx)){
    maf=exp(sqrt(mafR[idx]))
    tt=round((2^(seglogr.adj[idx]+1)-2*(1-rho))/rho,0)
    tt[tt<0]=0 #homdel
    mm=round((tt*maf*rho+(maf-1)*(1-rho))/(rho*(maf+1)),0)
    re=which(mm>tt) #rounding error can cause major>t
    if(any(re)){mm[re]=tt[re]}   
    t.em[idx]=tt
    major.em[idx]=mm

    minor.em[idx]=t.em[idx]-major.em[idx]
    genotype.em[idx] = paste("A",major.em[idx], "B", minor.em[idx], sep="")
    rhov.em[idx]=rho
    
    clonal.cluster[idx]=1
  }
  
  clonal.cluster[which.geno.em==4]=1
  
  
  #if het SNPs are too few, not sufficient information to estimate minor cn
  lownhet=which(nhet<min.nhet)
  minor.em[lownhet]=NA
  minor.em[t.em<=1]=0
  

  #set cf=1 for 2-1 segments (100% nothing)
  rhov.em[t.em==2&minor.em==1]=1
  rhov.em[t.em==2&is.na(minor.em)]=NA
  clonal.cluster[t.em==2&is.na(minor.em)]=NA
  clonal.cluster[chr==nX]=NA
  
  #for male, use the empirical call
  if(sum(chr==nX)>0){
    prop.nhet.chrX=sum(nhet[chr==nX])/sum(nmark[chr==nX])
    male=(prop.nhet.chrX<0.01)
  }else{
    male=FALSE
  }
  
  #normal male X is one copy. No het snps to start with, so don't call minor cn
  if(male){
    t.em[chr>=nX]=round(t.em[chr>=nX]/2,0)
    minor.em[chr>=nX]=NA
    normalX=which(t.em[chr>=nX]==1)
    if(any(normalX))rhov.em[chr>=nX][normalX]=1
  }
  
  
  
  out1=data.frame(seg[,1:9],start=startseq,end=endseq, cf.em=rhov.em,tcn.em=t.em, lcn.em=minor.em, clonal.cluster=clonal.cluster)

  if(rho<0.3){emflags=paste(emflags,"Low purity. Calls can be unreliable.",sep=" ")}

  out=list(purity=rho,ploidy=gamma,dipLogR=dipLogR,cncf=out1, emflags=emflags)
  
  return(out)
  
}

#####onepass###
onepass=function(x, trace, unif, rho, rhov, prior, posterior, sigma, min.nhet, rho.clust, maxiter, eps){  

  # added to make R CMD check stop complaining - didn't work
  # globalVariables(c("chr", "cnlr.median.clust", "dipLogR", "het", "jointseg", "logOR", "logOR2var", "logORvar", "logORvar.clust", "logR.adj", "mafR.clust", "major", "minor", "nX", "nclust", "ng", "nhet", "nmark", "nmark.clust", "segclust"))

  dif=1
  iter=0      
  
  while(dif>eps && iter<maxiter) {
    
    iter = iter + 1
    
    rhov[is.na(rhov)]=rho    
    #constraint: any segment cannot have purity higher than the mode purity 
    rhov[rhov>rho]=rho
    #avoid cf below 10%
    rhov[rhov<0.1]=rho
    
    rho.old = rho  
    rhov.old=rhov
    sigma.old = sigma
    prior.old=prior
    posterior.old=posterior
    
    ########
    #E-step#
    ########
    
    ####LogR mixture model parameter####
    gamma=2
    phi=2*(1-rho)+gamma*rho
    mu=log2(2*(1-rhov)+matrix(rhov,ncol=1)%*%t)-log2(phi)
    
    ####LogOR mixture model parameter####
    #allelic ratio
    k=(matrix(rhov,ncol=1)%*%major+1-rhov)/(matrix(rhov,ncol=1)%*%minor+1-rhov)
    logk=log(k)
    logk2=logk^2
    
    #posterior probability matrix
    #pmatrix=NULL
    pmatrix=matrix(NA,nrow=nrow(jointseg),ncol=ng)
    loglik=0
    
    clust=rep(segclust,nmark)
    segc=sort(unique(segclust[chr<=nX]))
    for(s in segc){
      idx=which(clust==s)
      x1ij=logR.adj[idx]
      upper=quantile(x1ij,0.95)
      lower=quantile(x1ij,0.05)
      x1ij[x1ij>upper]=NA
      x1ij[x1ij<lower]=NA
      mus=rep(mu[s,],each=length(idx))
      sd=sigma[s]
      if(rhov[s]<0.5){
        x1ij=rep(cnlr.median.clust[s]-dipLogR,length(idx))
        sd=0.1
      }
      #density for logR.adj (centered logR)
      d1=dnorm(x1ij,mean=mus,sd=sd)
      d1[d1==Inf]=NA
      
      #density for logOR, non-central chi-square
      nu=rep(logk2[s,],each=length(idx))
      lambda=nu/rep(logORvar[idx],ng)
      x2ij=logOR2var[idx]
      if(rhov[s]<0.5){
        x2ij=rep(mafR.clust[s]/logORvar.clust[s],length(idx))
        lambda=nu/logORvar.clust[s]
      }
      #d2=dchisq(x2ij+1,df=1,ncp=lambda)
      d2=dchisq(x2ij,df=1,ncp=lambda)
      d2=1/(abs(x2ij-lambda)+1e-6)
      d2[d2==Inf]=NA
      
      #likelihood
      d=d1*d2
      hetsum=d[rep(het[idx]==1,ng)]
      homsum=d1[rep(het[idx]==0,ng)]
      d=sum(hetsum[hetsum<Inf],na.rm=T)+sum(homsum[homsum<Inf],na.rm=T)
      if(!is.na(d)&d>0&d<Inf){loglik=loglik+log(d)}
      
      #heterozygous positions contribute to logR and logOR
      numerator1=matrix(d1*d2,nrow=length(idx),ncol=ng,byrow=F)
      numerator1=sweep(numerator1,MARGIN=2,prior[s,],`*`)
      
      #homozygous positions contribute to logR only
      numerator0=matrix(d1,nrow=length(idx),ncol=ng,byrow=F)
      numerator0=sweep(numerator0,MARGIN=2,prior[s,],`*`)
      
      numerator=numerator1
      numerator[het[idx]==0,]=numerator0[het[idx]==0,]
      
      tmp=apply(numerator,1,function(x)x/(sum(x,na.rm=T)+1e-5))
      #pmatrix=rbind(pmatrix,t(tmp))
      pmatrix[idx,]=t(tmp)
      
      #update prior
      prior[s,]=apply(t(tmp),2,function(x)mean(x,na.rm=T))
    }

    ########
    #M-step#
    ########
    
    #get CF per segments, pick mode close to 1 (favor high purity low cn solution)
    rhom=gammam=matrix(NA,nrow=nclust,ncol=ng)
    geno=matrix(0,nrow=nclust,ncol=ng)
    which.geno=posterior=rep(NA,nclust)
    for(i in segc){
      
      idx=which(clust==i)
      idxhet=which(clust==i&het==1)
      sump=apply(pmatrix[idx,,drop=F],2,function(x)sum(x,na.rm=T))
      
      #if probability is too small (highly uncertain), use lsd estimates for stability 
      if(all(is.na(prior[i,]))){
        prior[i,]=prior.old[i,]
      }else{
        if(sum(prior[i,],na.rm=T)==0)prior[i,]=prior.old[i,]
      }
      
      
      sump[prior[i,]<max(prior[i,])]=NA
      
      ##update k
      tmphet=pmatrix[idxhet,,drop=F]
      v1=as.vector((logOR[idxhet]^2-logORvar[idxhet])/logORvar[idxhet])
      v2=as.vector(1/logORvar[idxhet])
      sumdphet=apply(sweep(tmphet,MARGIN=1, v1, `*`), 2,function(x)sum(x,na.rm=T))
      sumphet=apply(sweep(tmphet,MARGIN=1,v2,`*`), 2,function(x)sum(x,na.rm=T))
      sumphet[is.na(sump)]=NA
      
      #CF from logOR    
      logk2hat=pmax(0,sumdphet/sumphet) #can be negative when k=1 logk=0 set to 0
      khat=exp(sqrt(logk2hat))
      a=(1-khat)/(khat*(minor-1)-(major-1))
      a[abs(a)==Inf]=NA
      a[a<=0]=NA
      a[a>1]=1
      if(all(nhet[segclust==i]<min.nhet))a=rep(NA,ng)
      
      #CF from logR
      tmp=pmatrix[idx,,drop=F]
      v=as.vector(logR.adj[idx])
      sumdp=apply(sweep(tmp,MARGIN=1,v,`*`),2,function(x)sum(x,na.rm=T))   
      mu.hat=sumdp/sump #mu.hat
      aa=2*(2^mu.hat-1)/(t-2)
      aa[abs(aa)==Inf]=NA
      aa[aa<=0]=NA
      aa[aa>1]=1
      
      aaa=pmin(a,aa,na.rm=T)
      #degenerate cases
      #homozygous deletion (0) and balanced gain (AABB, AAABBB), maf=0.5, purity information comes from logr only
      aaa[c(1,8,13)]=aa[c(1,8,13)]  
      
      #set upper bound at sample rho
      aaa=pmin(aaa,rho)
      
      #uniparental disomy (AA) CF information comes from logOR only.
      
      which.geno[i]=which.max(prior[i,])
      
      postprob=pmatrix[idx,which.geno[i]]
      posterior[i]=mean(postprob[postprob>0],na.rm=T)
      
      #update sigma
      y=as.vector(logR.adj[idx])*pmatrix[idx,,drop=F]
      r=y-mu[i,]*pmatrix[idx,,drop=F]
      ss=sqrt(sum(r[,which.geno[i]]^2,na.rm=T)/sum(pmatrix[idx,which.geno[i]])) 
      sigma[i]=ifelse(is.na(ss),0.5,ss)
      
      aaa[setdiff(1:ng,which.geno[i])]=NA  
      
      rhom[i,]=aaa 
      
    } #max prior
    
    #}

    #for segments noninformative for purity, plug in sample purity
    rhov=unlist(apply(rhom,1,function(x)ifelse(all(is.na(x)),NA,na.omit(x))))
    segclust.cf=rhov

    rhos=by(rep(rhov,as.vector(nmark.clust)),rep(rho.clust,as.vector(nmark.clust)),function(x)mean(x,na.rm=T))
    rhov=rhos[rho.clust]

    if(all(rho.clust==1))rhov=rep(rho,nclust)
    rhov=as.vector(rhov)
    rhov[rho.clust==1]=rho
    
    rho=max(rhov,na.rm=T)
    
    
    dif=c(abs(rhov-rhov.old),diff(posterior-posterior.old))
    dif = max(dif,na.rm=T)
    
    if(trace) {
      cat('iter:', iter, '\n')
      cat('dif:', dif, '\n')
      cat('purity:', rho, '\n')
      cat('ploidy:', gamma, '\n')
    }
    
  }
  
  out=list(posterior=posterior, which.geno=which.geno, rho=rho, rhov=rhov, segclust.cf=segclust.cf, rho.clust=rho.clust)
  
  return(out)
  
}
