#last update:07/21/15
#input
#x: output from procSample
#parameters
#mu: logR normal mean
#k: allelic ratio
#rho: tumor sample purity
#rhov: vector of cellular freuency of the alteration for each segment
#pmatrix, posterior probability matrix of dimension nseg (number of segments) by ng (number of genotypes)

emcncf=function(x,trace=FALSE,unif=FALSE,maxiter=10,eps=1e-3){  
  
  jointseg=x$jointseg
  out=x$out
  dipLogR=x$dipLogR
  jointseg=subset(jointseg,!is.na(cnlr))    
  logR=jointseg$cnlr 
  #segid=jointseg$segs
  logOR=jointseg$valor 
  logORvar=jointseg$lorvar  
  logOR2var=logOR^2/logORvar
  
  het=jointseg$het
  
  seg=out  
  nmark=seg$num.mark
  segs=rep(1:length(nmark),nmark)  
  nseg=length(nmark)
  if(nseg>300)stop("Likely hyper-segmented. Increase cval in procSample.")
  
  
  endseq=jointseg[cumsum(nmark),2]
  startseq=jointseg[c(1,cumsum(nmark)[-nrow(seg)]+1),2]
  seglen=(endseq-startseq)/1e6
  
  mafR=seg$mafR
  mafR[mafR<0]=0
  seglogr=seg$cnlr.median
  nclust=max(seg$segclust)
  
  emflags=NULL
  
  var=var(jointseg$cnlr,na.rm=T)    
  if(var>0.6){
    logR=rep(seglogr,nmark)
    emflags=paste(emflags,"Noisy sample, Calls can be unreliable.",sep=" ")
    #warning("Noisy sample. Calls unreliable.",call.=F)
  }
  
  #consider genotypes up to t=6, assume minor alelle is B, switch to cncf for high copy numbers (t>6) for computational efficiency
  genotype=c("0","A","AA","AB","AAB","AAA","AAAB","AABB","AAAA","AAAAB","AAABB","AAAAA","AAABBB","AAAABB","AAAAAB","AAAAAA")
  minor=c(0,0,0,1,1,0,1,2,0,1,2,0,3,2,1,0)
  major=c(0,1,2,1,2,3,3,2,4,4,3,5,3,4,5,6)
  t=ifelse(genotype=="0",0,nchar(genotype))
  
  ng=length(genotype)
  n=length(logR)  
    
  #if no big segment that is imbalanced then outplot diploid genome with purity=1
  if(max(mafR[seglen > 10], na.rm = T) < 0.05 & max(abs(seglogr)[seglen > 10], na.rm = T)<0.1){
    rhov.em=rep(1,nseg)
    major.em=rep(1,nseg); minor.em=rep(1,nseg)
    rho=NA
    gamma=2
    out1=data.frame(seg,cf.em=rhov.em,tcn.em=major.em+minor.em, lcn.em=minor.em)
    out1$tcn.em[is.na(out1$tcn.em)] <- out1$tcn[is.na(out1$tcn.em)]
    emflags=paste(emflags,"Insufficient information. Likely diplod or purity too low.",sep=" ")
    out=list(purity=rho,ploidy=gamma,cncf=out1,emflags=emflags)
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
  genotype.lsd[homdel]="0"
  
  a=lapply(1:sum(!nas),function(x)paste(rep("A",major.lsd[!nas][x]),collapse=""))
  b=lapply(1:sum(!nas),function(x)paste(rep("B",minor.lsd[!nas][x]),collapse=""))
  genotype.lsd[!nas]=unlist(lapply(1:sum(!nas),function(x)paste(a[[x]],b[[x]],sep="")))
  which.geno.lsd=match(genotype.lsd,genotype)
    
  rhov.lsd[major.lsd==1&minor.lsd==1]=NA
  #rhov.lsd[nas]=NA
  meanrhov.lsd=mean(rhov.lsd,na.rm=T)
  rhov.lsd.subset=rhov.lsd
  rhov.lsd.subset[which.geno.lsd%in%c(5,7,10,11,14,15)]=NA 
  rhov.lsd.subset[t.lsd>6]=NA
  #rhov.lsd[(major.lsd!=minor.lsd)&minor.lsd!=0]=NA #imbalanced seg has big identifiability issue
  loh=which(major.lsd %in% c(1,2)& minor.lsd==0 & seglen>10) #use only LOH seg for initial estimate
  if(length(loh)>1){
    rho=max(rhov.lsd[loh],na.rm=T)
  }else{  
    rho=mean(rhov.lsd.subset,na.rm=T)
  }
  if(is.na(rho)){rho=meanrhov.lsd}
  
  lowpur=FALSE  
  #if(rho<0.35){lowpur=TRUE}
  
  if(lowpur){
    #rho=max(0.5,rho)
    rhov=rep(rho,nclust)
    }else{
  rhov=rhov.lsd   
  rhov[is.na(rhov)]=rho
  #avoid initial value too low
  #rhov[rhov<0.1]=rho
  #avoid 1
  rhov[rhov==1]=rho  
  rhov=as.vector(by(rhov,seg$segclust,mean))
  }

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
    logk2=log(k)^2
    
    #posterior probability matrix
    pmatrix=NULL
    loglik=0
    
    clust=rep(seg$segclust,seg$num.mark)
    for(s in 1:nclust){
      idx=which(clust==s)
      
      #density for logR.adj (centered logR)
      d1=dnorm(logR.adj[idx],mean=rep(mu[s,],each=length(idx)),sd=sigma[s])
      d1[d1==Inf]=NA
      
      #density for logOR, non-central chi-square
      lambda=rep(logk2[s,],each=length(idx))
      lambda=lambda/rep(logORvar[idx],ng)
      d2=dchisq(logOR2var[idx],df=1,ncp=lambda)
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
      
      tmp=apply(numerator,1,function(x)x/(sum(x)+1e-5))
      pmatrix=rbind(pmatrix,t(tmp))
      
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
    for(i in 1:nclust){
      
      idx=which(clust==i)
      idxhet=which(clust==i&het==1)
      sump=apply(pmatrix[idx,,drop=F],2,function(x)sum(x,na.rm=T))
      
      #if probability is too small (highly uncertain), use lsd estimates for stability 
      if(all(is.na(prior[i,]))){prior[i,]=prior.old[i,]}
      if(max(prior[i,],na.rm=T)>0.05){   
        
      #calculate rho for the most likely genotype(s) for segment i
      #if there more more than one likely candidates save two and pick one with higher CF
      top2=sort(prior[i,],decreasing=T)[1:2]
      if(top2[2]>0.05&abs(diff(top2))<0.001){
          sump[prior[i,]<quantile(prior[i,],(ng-2)/ng)]=NA
       }else{
          sump[prior[i,]<max(prior[i,])]=NA
       }
        
        ##update k
        tmphet=pmatrix[idxhet,,drop=F]
        sumdphet=apply(sweep(tmphet,MARGIN=1,as.vector((logOR[idxhet]^2-logORvar[idxhet])/logORvar[idxhet]),`*`), 2,function(x)sum(x,na.rm=T))
        sumphet=apply(sweep(tmphet,MARGIN=1,as.vector(1/logORvar[idxhet]),`*`), 2,function(x)sum(x,na.rm=T))
        sumphet[is.na(sump)]=NA
                
        #CF from logOR    
        logk2hat=pmax(0,sumdphet/sumphet) #can be negative when k=1 logk=0 set to 0
        khat=exp(sqrt(logk2hat))
        a=(1-khat)/(khat*(minor-1)-(major-1))
        a[abs(a)==Inf]=NA
        a[a<=0|a>=1]=NA
        
        #CF from logR
        tmp=pmatrix[idx,,drop=F]
        sumdp=apply(sweep(tmp,MARGIN=1,as.vector(logR.adj[idx]),`*`),2,function(x)sum(x,na.rm=T))   
        mu.hat=sumdp/sump #mu.hat
        aa=2*(2^mu.hat-1)/(t-2)
        aa[abs(aa)==Inf]=NA
        aa[aa<=0|aa>=1]=NA
        
        aaa=pmax(a,aa,na.rm=T)
        #degenerate cases
        #homozygous deletion (0) and balanced gain (AABB, AAABBB), maf=0.5, purity information comes from logr only
        aaa[c(1,8,13)]=aa[c(1,8,13)]  
     
        #set upper bound at sample rho
        aaa=pmin(aaa,rho)
        
        #het dip (AB) seg has no information, set CF=0.8 
        if(any(which(!is.na(sump))==4)){aaa[4]=0.8}
        
        #uniparental disomy (AA) CF information comes from logOR only.
        
        #if there are two likely genotype, choose one with higher purity (e.g.,AAB 80% or AAAB 50%)
        #if the higher CF exceeds sample purity, then the lower CF is the right one
        if(all(is.na(aaa))){which.geno[i]=which.max(prior[i,])}else{
          which.geno[i]=ifelse(max(aaa,na.rm=T)<rho,which.max(aaa),which.min(aaa))
        }
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
      
    }
    
    
    #for segments noninformative for purity, plug in sample purity
    rhov=unlist(apply(rhom,1,function(x)ifelse(all(is.na(x)),NA,na.omit(x))))

    #Determine sample rho
    #pick mode closest to 1
    rhov.long=rhov[seg$segclust]
    which.geno.long=which.geno[seg$segclust] 
    rhov.long[which.geno.long == 4]=NA
    #EM tend to overfit for homdel
    if(sum(which.geno.long == 1,na.rm=T)>0){
      idx=which(which.geno.long == 1)
      rhov.long[idx]=rhov.lsd[idx]
      which.geno.long[idx]=which.geno.lsd[idx]
    }
    meanrho=mean(rhov.long,na.rm=T)
    rhov.long.subset=rhov.long
    rhov.long.subset[which.geno.long %in% c(5,7,10,11,14,15)]=NA #Imbalanced gains have big identifiability issue
    
    #if low purity, assume rhov the same
    #if more than 100 seg give rho estimate, use density to find rho
    #otherwise use LOH seg for purity estimate
    if(all(is.na(rhov.long))){
      rho=max(rhov.lsd,na.rm=T)
      }else{  
       if(lowpur){          
          rho=mean(rhov.long[rhov.long>0.1],na.rm=T)     
          rhov=rep(rho,nclust)      
        }else{       
            nona=na.omit(rhov.long)        
            if(length(nona)>100){
               rho=find.mode(nona)$rho
            }else{
              loh=rhov.long[major[which.geno.long]%in%c(1,2) & minor[which.geno.long]==0 & seglen>10]
              if(length(loh)==0){rho=mean(rhov.lsd.subset,na.rm=T)}else{
               rho=max(loh,na.rm=T)
              }
           }
      if(is.na(rho)){rho=meanrho}        
      #for noninformative region, plug in sample purity
      rhov[is.na(rhov)]=rho    
      #constraint: any segment cannot have purity higher than the mode purity 
      rhov[rhov>rho]=rho
      #rhov[rhov<0.1]=rho
        }
    }
  
    
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
  
  #no information from A1B1 segments
  rhov.em[which.geno.em==4]=NA
  
  
  genotype.em=which.geno.em  
  genotype.em[which(!is.na(which.geno.em))]=paste("A",major[which.geno.em],"B",minor[which.geno.em],sep="")[which(!is.na(which.geno.em))]
  t.em=t[which.geno.em]
  major.em=major[which.geno.em]
  minor.em=minor[which.geno.em]
  
  #calculate ploidy
  gamma=(2^(-dipLogR)*(2*(1-rho)+2*rho)-2*(1-rho))/rho
  
  #hybrid: for high copy number (t>6), swicth to lsd estimate
  seglogr.adj=seg$cnlr.median-dipLogR
  idx=which(seglogr.adj>1.6*rho|is.na(which.geno.em))
  if(any(idx)){
    maf=exp(sqrt(seg$mafR[idx]))
    t.em[idx]=round((2^(seglogr.adj[idx]+1)-2*(1-rho))/rho,0)
    major.em[idx]=round((t.em[idx]*maf*rho+(maf-1)*(1-rho))/(rho*(maf+1)),0)
    minor.em[idx]=t.em[idx]-major.em[idx]
    genotype.em[idx] = paste("A",major.em[idx], "B", minor.em[idx], sep="")
    rhov.em[idx]=rho
    
    #genotype.em[idx]=genotype.lsd[idx]
    #t.em[idx]=t.lsd[idx]
    #minor.em[idx]=minor.lsd[idx]
    #major.em[idx]=major.lsd[idx]

  }
  
  #EM over-calling homozygous deletion, switch to lsa
  idx=which(which.geno.em==1&seglen>10)
  if(any(idx)){
    genotype.em[idx] = "AB"
    t.em[idx]=2
    minor.em[idx]=1
    major.em[idx]=1
    rhov.em[idx]=NA    
    
    #genotype.em[idx]=genotype.lsd[idx]
    #t.em[idx]=t.lsd[idx]
    #minor.em[idx]=minor.lsd[idx]
    #major.em[idx]=major.lsd[idx]
  }

  out1=data.frame(seg,cf.em=rhov.em,tcn.em=major.em+minor.em, lcn.em=minor.em)
  out1$tcn.em[is.na(out1$tcn.em)] <- out1$tcn[is.na(out1$tcn.em)]
  
  if(rho<0.3){emflags=paste(emflags,"Low purity. Calls can be unreliable.",sep=" ")}

  out=list(loglik=loglik,purity=rho,ploidy=gamma,dipLogR=dipLogR, seglen=seglen,cncf=out1, emflags=emflags)
  
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

plotSampleCNCF=function(x,fit) {
  mat=x$jointseg
  mat=subset(mat,chrom<23)
  
  out=subset(x$out,chrom<23)
  
  cncf=fit$cncf
  cncf=subset(cncf,chrom<23)
     
  dipLogR <- fit$dipLogR
  
  layout(matrix(c(1,1,2,2,3,3,4,4,5,5,6,6), ncol=1))
  par(mar=c(3,3,1,1), mgp=c(2, 0.7, 0))

  chr=mat$chrom
  len=table(chr)
  altcol=rep(c("light blue","gray"),12)[-c(23:24)]
  chr.col=rep(altcol,len)
  nmark=cncf$num.mark
  tmp=cumsum(len)
  start=c(1,tmp[-22]+1)
  end=tmp
  mid=start+len/2
  
  
  plot(mat$cnlr, pch=".", axes=F,cex=1.5,ylim=c(-3,3), col = c("grey","lightblue")[1+rep(cncf$chrom-2*floor(cncf$chrom/2), cncf$num.mark)], ylab="log-ratio")
  abline(h=dipLogR, col="magenta4")
  points(rep(cncf$cnlr.median, cncf$num.mark), pch=".", cex=2, col="brown")
  axis(side=1,at=mid,1:22,cex.axis=1,las=2)
  axis(side=2)
  box()
  
  
  plot(mat$valor, axes=F, pch=".", cex=1.5, col = c("grey","lightblue")[1+rep(cncf$chrom-2*floor(cncf$chrom/2), cncf$num.mark)], ylab="log-odds-ratio", ylim=c(-4,4))
  points(rep(sqrt(abs(cncf$mafR)), cncf$num.mark), pch=".", cex=2, col="brown")
  points(-rep(sqrt(abs(cncf$mafR)), cncf$num.mark), pch=".", cex=2, col="brown")
  axis(side=1,at=mid,1:22,cex.axis=1,las=2)
  axis(side=2)
  box()
    
  plot(rep(cncf$cf.em, cncf$num.mark), axes=F,pch=".", cex=2, xlab="Chromosome",ylab="Cellular fraction (EM)", ylim=c(0,1))
  axis(side=1,at=mid,1:22,cex.axis=1,las=2)
  axis(side=2,cex=0.8)
  box()
  abline(v=start,lty=3,col="gray")
  abline(v=end,lty=3,col="gray")
  #abline(h=c(0.2,0.4,0.6,0.8),lty=3,col="gray")
    
  
  # scale tcn so that very high copy numbers don't take up space
  tcnscaled <- cncf$tcn.em
  tcnscaled[cncf$tcn.em > 5 & !is.na(cncf$tcn.em)] = (5 + (tcnscaled[cncf$tcn.em > 5& !is.na(cncf$tcn.em)] -5)/3)
  matplot(cbind(rep(tcnscaled, cncf$num.mark), rep(cncf$lcn.em,cncf$num.mark)-0.1), pch=".", cex=3, col=1:2, lwd=1, ylab="Integer copy number (EM)", yaxt="n",xaxt="n")
  axis(2, at=c(0:5,5+(1:35)/3), labels=0:40)
  axis(side=1,at=mid,1:22,cex.axis=1,las=2)
  box()
  abline(v=start,lty=3,col="gray")
  abline(v=end,lty=3,col="gray")
  abline(h=c(0:5,5+(1:35)/3),lty=3,col="gray")
  
  
  plot(rep(cncf$cf, cncf$num.mark), axes=F,pch=".", cex=2, xlab="Chromosome",ylab="Cellular fraction (cncf)", ylim=c(0,1))
  axis(side=1,at=mid,1:22,cex.axis=1,las=2)
  axis(side=2,cex=0.8)
  box()
  abline(v=start,lty=3,col="gray")
  abline(v=end,lty=3,col="gray")
  #abline(h=c(0.2,0.4,0.6,0.8),lty=3,col="gray")
  
  tcnscaled <- cncf$tcn
  tcnscaled[cncf$tcn > 5 & !is.na(cncf$tcn)] = (5 + (tcnscaled[cncf$tcn > 5& !is.na(cncf$tcn)] -5)/3)
  matplot(cbind(rep(tcnscaled, cncf$num.mark), rep(cncf$lcn,cncf$num.mark)-0.1), pch=".", cex=3, col=1:2, lwd=1, ylab="Integer copy number (cncf)", yaxt="n",xaxt="n")
  axis(2, at=c(0:5,5+(1:35)/3), labels=0:40)
  axis(side=1,at=mid,1:22,cex.axis=1,las=2)
  box()
  abline(v=start,lty=3,col="gray")
  abline(v=end,lty=3,col="gray")
  abline(h=c(0:5,5+(1:35)/3),lty=3,col="gray")
    
}


