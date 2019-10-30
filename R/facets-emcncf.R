#genotype mixture model using EM algorithm to call allele-specific copy number and cellular fraction
emcncf=function(x,trace=FALSE,unif=FALSE,min.nhet=15,maxiter=10,eps=1e-3){  
  
  jointseg=x$jointseg
  out=x$out
  dipLogR=x$dipLogR
  nX=x$nX
  seg=out
  
  jointseg=subset(jointseg,!is.na(jointseg$cnlr))  
  logR=jointseg$cnlr 
  logOR=jointseg$valor 
  logORvar=jointseg$lorvar 
  logOR2var=logOR^2/logORvar
  logORvar.clust=by(logORvar,jointseg$segclust,function(x)mean(na.omit(x)))
  
  het=jointseg$het
  
  nmark=seg$num.mark
  segclust=seg$segclust
  cnlr.median.clust=by(seg$cnlr.median,segclust,function(x)mean(na.omit(x)))
  #mafR.clust=by(seg$mafR,segclust,function(x)mean(na.omit(x)))
  mafR.clust=by(seg$mafR.clust,segclust,function(x)mean(na.omit(x)))
  segs=rep(1:length(nmark),nmark)  
  nseg=length(nmark)
  nhet=seg$nhet
  chr=seg$chrom
  #if(nseg>500)stop("Likely hyper-segmented. Increase cval in procSample.")
  
  
  endseq=jointseg[cumsum(nmark),2]
  startseq=jointseg[c(1,cumsum(nmark)[-nrow(seg)]+1),2]
  seglen=(endseq-startseq)/1e6
  
  mafR=seg$mafR
  mafR[mafR<0]=0
  seglogr=seg$cnlr.median
  nclust=max(segclust)
  
  emflags=NULL
  
  var=var(jointseg$cnlr,na.rm=T)    
  if(var>0.6){
    logR=rep(seglogr,nmark)
    emflags=paste(emflags,"Noisy sample, Calls can be unreliable.",sep=" ")
  }
  
  #consider genotypes up to t=6, assume minor alelle is B, switch to cncf for high copy numbers (t>6) for computational efficiency
  genotype=c("0","A","AA","AB","AAB","AAA","AAAB","AABB","AAAA","AAAAB","AAABB","AAAAA","AAABBB","AAAABB","AAAAAB","AAAAAA")
  minor=c(0,0,0,1,1,0,1,2,0,1,2,0,3,2,1,0)
  major=c(0,1,2,1,2,3,3,2,4,4,3,5,3,4,5,6)
  t=ifelse(genotype=="0",0,nchar(genotype))
  
  ng=length(genotype)
  n=length(logR)  
    
  #diploid genome with purity=1
  if(all(seg$tcn[seg$chrom<nX]==2 & seg$lcn[seg$chrom<nX]%in%c(1, NA))|max(mafR.clust[seg$chrom<nX & seg$nhet>min.nhet], na.rm = T) < 0.05){
    rho=NA
    gamma=2
    out1=data.frame(seg[,1:9],start=startseq,end=endseq,cf.em=seg$cf,tcn.em=seg$tcn,lcn.em=seg$lcn)
    emflags=paste(emflags,"Insufficient information to estimate purity. Likely diplod or purity too low.",sep=" ")
    out=list(purity=rho,ploidy=gamma,dipLogR=dipLogR,cncf=out1, emflags=emflags)
    return(out)    
    stop("Insufficient information",call.=F)
  }

  
  #intialize cellular fraction rho vector use least squared distance estimates
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
  
  #if(!all(is.na(rhov.lsd[seglen>35]))){
  #  naive=max(by(rhov.lsd[seglen>35],segclust[seglen>35],function(x)mean(x,na.rm=T)),na.rm=T)
  #}else{
  naive=quantile(rhov.lsd,prob=0.75,na.rm=T)
  #}

  
  rhov.lsd.subset=rhov.lsd
  rhov.lsd.subset[which.geno.lsd%in%c(3,5,7,10,11,14,15,NA)]=NA 
  rhov.lsd.subset[t.lsd>6]=NA
  rhov.lsd.subset[seglen<50]=NA
  loh=which(t.lsd>=1 & minor.lsd==0 & seglen>50) #use only LOH seg for initial estimate
  
  rho=NA
  if(length(loh)>2){
    rho=max(by(rhov.lsd[loh],segclust[loh],function(x)mean(x,na.rm=T)),na.rm=T)
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
  rhov=as.vector(by(rhov,segclust,mean))
  rhov0=rhov
  
  lowpur=FALSE
  #if(rho<0.2){lowpur=TRUE}
  #if(lowpur){
  #rhov=rep(rho,nclust)
  #}

  #center logR at diphet
  logR.adj=logR-dipLogR
    
  #initial value for genotype priror
  prior=matrix(1/ng,nrow=nclust,ncol=ng)
  
  #initial value for sigma 
  sigma=rep(2,nclust)  
  
  #cold start for rho if set unif=TRUE
  if(unif){
    rhov=runif(nclust,0.3,0.8)
  }
  
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
      if(rhov[s]<0.4){
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
      if(rhov[s]<0.4){
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
      
      if(max(prior[i,],na.rm=T)>0.05){   
        
      ##calculate rho for the most likely genotype(s) for segment i
      ##if there more more than one likely candidates save two and pick one with higher CF
      #top2=sort(prior[i,],decreasing=T)[1:2]
      #if(top2[2]>0.05&abs(diff(top2))<0.0001){
      #    sump[prior[i,]<quantile(prior[i,],(ng-2)/ng)]=NA
      # }else{
      #    sump[prior[i,]<max(prior[i,])]=NA
      # }
      
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
        
        aaa=pmax(a,aa,na.rm=T)
        #degenerate cases
        #homozygous deletion (0) and balanced gain (AABB, AAABBB), maf=0.5, purity information comes from logr only
        aaa[c(1,8,13)]=aa[c(1,8,13)]  
     
        #set upper bound at sample rho
        aaa=pmin(aaa,rho)
        
        #uniparental disomy (AA) CF information comes from logOR only.
        
        ##if there are two likely genotype, choose one with higher purity (e.g.,AAB 80% or AAAB 50%)
        ##if the higher CF exceeds sample purity, then the lower CF is the right one
        #if(all(is.na(aaa))){which.geno[i]=which.max(prior[i,])}else{
        #  which.geno[i]=ifelse(max(aaa,na.rm=T)<rho,which.max(aaa),which.min(aaa))
        #}
        
        which.geno[i]=which.max(prior[i,])
        
        postprob=pmatrix[idx,which.geno[i]]
        posterior[i]=mean(postprob[postprob>0],na.rm=T)
        
        #update sigma
        y=as.vector(logR.adj[idx])*pmatrix[idx,,drop=F]
        r=y-mu[i,]*pmatrix[idx,,drop=F]
        ss=sqrt(sum(r[,which.geno[i]]^2,na.rm=T)/sum(pmatrix[idx,which.geno[i]])) 
        sigma[i]=ifelse(is.na(ss),0.5,ss)
        
        aaa[setdiff(1:ng,which.geno[i])]=NA  
        
        #het dip (AB) seg has no information, set CF at a high value less than 1
        #if(any(which(!is.na(sump))==4)){aaa[4]=0.9}
        
        rhom[i,]=aaa 
        
      } #max prior
      
    }
    
    #for segments noninformative for purity, plug in sample purity
    rhov=unlist(apply(rhom,1,function(x)ifelse(all(is.na(x)),NA,na.omit(x))))

    #Determine sample rho
    rhov.long=rhov[seg$segclust]
    which.geno.long=which.geno[seg$segclust] 
    rhov.long[which.geno.long == 4]=NA
    rhov.long[chr>=nX]=NA
    
    if(sum(!is.na(rhov.long[seglen>35]))>1){
     meanrho=max(by(rhov.long[seglen>35],segclust[seglen>35],function(x)mean(x,na.rm=T)),na.rm=T)
    }else{
     meanrho=quantile(rhov.long,prob=0.75,na.rm=T) #if no big segments, probably very noisy sample or over-segemnted. can't estimate rho accurately, just take upper quatile
    }
    
    rhov.long.subset=rhov.long
    rhov.long.subset[which.geno.long %in% c(5,7,10,11,14,15,NA)]=NA #Imbalanced gains have big identifiability issue
    rhov.long.subset[seglen<50]=NA
    
    #if low purity, assume rhov the same
    #if more than 100 seg give rho estimate, use density to find rho
    #otherwise use LOH seg for purity estimate
    if(all(is.na(rhov.long))){
      rho=naive
      }else{  
       if(lowpur){          
          rho=max(rhov.long,na.rm=T)     
          rhov=rep(rho,nclust)      
        }else{       
            nona=na.omit(rhov.long)        
            if(length(nona)>100){
               rho=find.mode(nona)$rho
            }else{
              loh=which(major[which.geno.long]>=1 & minor[which.geno.long]==0 & seglen>50)
              rhov.loh=rep(NA,length(rhov.long))
              rhov.loh[loh]=rhov.long[loh]
              if(length(loh)>1 & !all(is.na(rhov.loh))){
                rho=max(by(rhov.loh,segclust,function(x)mean(x,na.rm=T)),na.rm=T)
              }else{
                if(length((na.omit(rhov.long.subset)))>1&dipLogR< -0.3)rho=max(by(rhov.long.subset,segclust,function(x)mean(x,na.rm=T)),na.rm=T)
              }
            }
          }
         }
  
    if(is.na(rho))rho=meanrho
    dif = quantile(abs(rhov-rhov.old),0.9,na.rm=T)

    if(trace) {
      cat('iter:', iter, '\n')
      cat('dif:', dif, '\n')
      cat('purity:', rho, '\n')
      cat('ploidy:', gamma, '\n')
    }
    
  }
  
  rhov.em=rhov[seg$segclust]
  which.geno.em=which.geno[seg$segclust]
  
  genotype.em=which.geno.em  
  genotype.em[which(!is.na(which.geno.em))]=paste("A",major[which.geno.em],"B",minor[which.geno.em],sep="")[which(!is.na(which.geno.em))]
  t.em=t[which.geno.em]
  major.em=major[which.geno.em]
  minor.em=minor[which.geno.em]
  
  #calculate ploidy
  gamma=(2^(-dipLogR)*(2*(1-rho)+2*rho)-2*(1-rho))/rho
  
  #hybrid: for high copy number (t>6), use moment estimates
  seglogr.adj=seg$cnlr.median-dipLogR
  idx=which(seglogr.adj>1.6*rho|is.na(which.geno.em))
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

  }
  

#  #EM over-calling homozygous deletion at low cf, switch to lsa
#   idx=which(which.geno.em==1&rhov.em<0.25&seglen>35)
#   if(any(idx)){
#     genotype.em[idx]=genotype.lsd[idx]
#     t.em[idx]=t.lsd[idx]
#     minor.em[idx]=minor.lsd[idx]
#     major.em[idx]=major.lsd[idx]
#     rhov.em[idx]=rhov.lsd[idx]
#   }
  
  
  #if het SNPs are too few, not sufficient information to estimate minor cn
  lownhet=which(nhet<min.nhet)
  minor.em[lownhet]=NA
  minor.em[t.em<=1]=0

  #set cf=1 for 2-1 segments (100% nothing)
  rhov.em[t.em==2&minor.em==1]=1
  rhov.em[t.em==2&is.na(minor.em)]=NA
  
  
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
  
  #out1=data.frame(seg,cf.em=rhov.em,tcn.em=t.em, lcn.em=minor.em)
  out1=data.frame(seg[,1:9],start=startseq,end=endseq,cf.em=rhov.em,tcn.em=t.em,lcn.em=minor.em)

  if(rho<0.3){emflags=paste(emflags,"Low purity. Calls can be unreliable.",sep=" ")}

  out=list(loglik=loglik,purity=rho,ploidy=gamma,dipLogR=dipLogR,seglen=seglen,cncf=out1, emflags=emflags)
  
  return(out)
  
}

find.mode=function (x) 
{
  den = density(na.omit(x),n=length(x))
  y = den$y
  y[y < 0.001] = 0
  rho = den$x[which.max(y)]
  difseq = rle(sign(diff(y)))
  nmodes = length(difseq$lengths)/2
  len = cumsum(difseq$lengths)
  modes.pos = len[which(difseq$values == 1)] + 1
  modes = y[modes.pos]
  idx = order(modes, decreasing = T)
  signodes = modes.pos[which(modes > 1.5)]
  if (length(signodes) >= 2) {
    #rho = max(den$x[modes.pos[idx]],na.rm=T)
    rho = max(den$x[signodes],na.rm=T)
  }
  out = list(rho = rho, nmodes = nmodes)
  out
}

# added emcncf2 to defunct
emcncf2 <- function() {
    .Defunct("emcncf", msg="emcncf2 is defunct as of v0.6.0")
}
