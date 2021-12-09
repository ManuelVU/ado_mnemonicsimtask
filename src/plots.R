# This code can be used to generate the figures in the "Adaptive Design 
# Optimization for a Mnemonic Similarity Task" Figures in the paper are 
# numbered in the code, however, it does not generate particular examples but 
# figures for all participants when needed. Figures 2 and 3 are examples and 
# are generated independently.

#### Load data and libraries needed ####
library(R.matlab)
parameters <- 7
times <- 192
trials <- 192
n.sub <- 20
ages <- 2

# load color pallete
colors_pantone <- readMat('src/PantoneFall2016.mat')
colors.p <- c()
for(i in 1:10){
  colors.p[i] <- rgb(colors_pantone$pantone[[i]], maxColorValue=1)
}

# assign 3 colors with transparency, colors are yellow (experimental designs), blue (young) and red (elderly)
col.a <- rev(c(paste(c(colors.p[4],'33'),collapse=''),
               paste(c(colors.p[1],'33'),collapse=''),
               paste(c(colors.p[8],'33'),collapse='')))

# Assign the same 3 colors without transparency
col.al <- rev(colors.p[c(4,1,8)])

#### Figure 1: Proportion of new responses by age group ####
# estimate the proportion of 'new' responses by participant in each population and for each stimulus type
load("data/memory.RData")
p.n <- array(NA,dim=c(n.sub,parameters,ages))
for(a in 1:ages){
  for(i in 1:n.sub){
    for(s in 1:parameters){
      p.n[i,s,a] <- sum(results$res[results$st[,i,a]==s,i,a]==1)/length(results$res[results$st[,i,a]==s,i,a]==1)
    }
  }  
}

pdf('figures/proportionofcorrect.pdf',width = 9,height = 8)
par(oma=c(4.5,4.5,0.1,0.1),
    mai=c(0.1,0.1,0.1,0.1))
centers <- rbind(seq(0,6)-0.125,seq(0,6)+0.125)
plot(0,0,xlim=c(-0.2,6.2),ylim=c(-0.001,1.1),type="n",axes=F,ann=F)
box(bty='l')
axis(1, at=seq(0,6),labels = c('criterion','lure 1','lure 2','lure 3','lure 4','lure 5','new'),
     cex.axis=1.5)
axis(2,las=2,cex.axis=1.5,at=seq(0,1,0.2),labels=c('0',seq(0.2,0.8,0.2),'1'))
for(a in 1:2){
  for(s in 1:7){
    qx <- quantile(p.n[,s,a],probs = c(0.025,0.25,0.5,0.75,0.975))
    rect(xleft = centers[a,s]-0.1  ,xright = centers[a,s]+0.1,ybottom = qx[2], 
         ytop = qx[4], col = col.a[(a+1)],
         border = col.al[(a+1)])
    segments(x0=rep(centers[a,s],2),x1=rep(centers[a,s],2),
             y0=c(qx[1],qx[4]),y1=c(qx[2],qx[5]),col=col.al[(a+1)],lwd=1.5)
    points(centers[a,s],qx[3],pch=16,col=col.al[(a+1)],cex=1.3)
    points(jitter(rep(centers[a,s],dim(p.n)[1]),amount = 0.07),p.n[,s,a],
           pch=16,col=col.a[(a+1)],cex=1.3)
  }
}
mtext('Proportion of New Responses',cex=1.7,side=2,line=3.3)
legend('topleft',legend = c("Young","Elderly"),pch=16,col=col.al[c(2,3)],bty='n',
       cex=1.4)
dev.off()

#### Figure 4 and Figure 5: k and d' by participant, examples used are young 15 (fig. 4) and elderly 14 (fig.5) ####
load('data/contaminant_sdt.RData')
ylable <- c('Criterion','Lure 1','Lure 2','Lure 3','Lure 4','Lure 5','New')
age.name <- c('Young','Elderly')
for(aa in 1:ages){
  for(pp in 1:n.sub){
    file.name <- paste(c('individualD_',bquote(.(age.name[aa])),'_',bquote(.(pp)),'.pdf'),collapse='')
    pdf(file=paste(c('figures/individualdprime/',file.name),collapse=''))
    par(oma=c(2,2,2,.7))
    layout(rbind(c(1,2),
                 c(3,4),
                 c(5,6),
                 c(7,8)))
    par(mai=c(0.3,0.5,0.2,0.1),xaxs="i",yaxs='i')
    counts <- c()
    for(ss in 1:7){
      counts <- append(counts,which(sdt$optimal$storder[,1,pp,aa]==ss)[-c(1:table(results$st[,pp,aa])[ss])])
      plot(0,0,xlim=c(-1,trials),type="n",axes=F,ann=F,
           ylim=c(min(sdt$data$ci[,1,,pp,aa],sdt$optimal$ci[,1,,pp,aa]),
                  max(sdt$data$ci[,2,,pp,aa],sdt$optimal$ci[,2,,pp,aa])))
      if(ss==7){
        mtext(expression(paste("d'")),side=2,las=2,line=2.3,cex=1.5)
        mtext("Trial",side=1,line=3,cex=1.5)
      }
      mtext(paste(ylable[ss]),side=3)
      axis(1,at=c(1,50,100,150,192),cex.axis=1.25)
      axis(2,las=2,cex.axis=1.25)
      polygon(x=c(seq(1,times),rev(seq(1,times))),border=F,col=col.a[1],
              y=c(sdt$data$ci[ss,2,,pp,aa],rev(sdt$data$ci[ss,1,,pp,aa])))
      lines(seq(1,times),sdt$optimal$ci[ss,1,,pp,aa],col=col.al[aa+1],lwd=1.2)
      lines(seq(1,times),sdt$optimal$ci[ss,2,,pp,aa],col=col.al[aa+1],lwd=1.2)
    }
    plot(0,0,ylim=c(min(sdt$data$ut[,pp,aa],sdt$optimal$ut[,pp,aa]),
                    max(sdt$data$ut[,pp,aa],sdt$optimal$ut[,pp,aa])),
         type="n",axes=F,ann=F,xlim=c(0,dim(sdt$data$ut)[1]))
    legend('topright',bty='n',col=c(col.al[1],col.al[aa+1]),
           legend=c('Experiment','ADO'),cex=1.3,
           pch=c(16,16))
    lines(smooth.spline(seq(1,dim(sdt$data$ut)[1]),sdt$data$ut[,pp,aa]),col=col.al[1],lwd=2)
    lines(smooth.spline(seq(1,dim(sdt$optimal$ut)[1]),sdt$optimal$ut[,pp,aa]),col=col.al[aa+1],lwd=2)
    points(seq(1,dim(sdt$data$ut)[1]),sdt$data$ut[,pp,aa],col=col.a[1],pch=16,cex=0.8)
    points(seq(1,dim(sdt$optimal$ut)[1]),sdt$optimal$ut[,pp,aa],col=col.a[aa+1],pch=16,cex=0.8)
    axis(1,at=c(1,50,100,150,192),cex.axis=1.25)
    axis(2,las=2,cex.axis=1.25)
    mtext("KL Divergence",side=2,line=2.7,cex=1.1)
    abline(v=min(counts),col='slategray',lty=2)
    dev.off()
  }
}

#### Figure 6: KL quantiles by age group ####
load('data/contaminant_sdt.RData')
klq <- array(NA,dim=c(3,192,2,2))
# obtain the quantiles of the KL divergence for each age group and each trial
klq[1:3,,1,1] <- apply(sdt$data$ut[,,1],1,quantile)[c(2,3,4),]
klq[1:3,,2,1] <- apply(sdt$data$ut[,,2],1,quantile)[c(2,3,4),]
klq[1:3,,1,2] <- apply(sdt$optimal$ut[,,1],1,quantile)[c(2,3,4),]
klq[1:3,,2,2] <- apply(sdt$optimal$ut[,,2],1,quantile)[c(2,3,4),]

pdf(file='figures/quantile_KL.pdf',width=8,height = 4)
par(oma=c(2.3,3,1,1))
layout(t(c(1,2)))
par(mai=c(0.3,0.3,0.2,0.1),xaxs="i",yaxs='i')
plot(seq(1,192),ylim=c(0,max(klq[,,1,])),type='n',axes=F,ann=F)
box(bty='l')
axis(1,at=c(1,50,100,150,192))
axis(2,las=2)
polygon(c(seq(1,192),rev(seq(1,192))),
        c(klq[3,,1,1],rev(klq[1,,1,1])),
        col=col.a[1],border=F)
mtext('KL Divergence',side=2,line=2.5,cex=1.3)
lines(seq(1,192),klq[2,,1,1],col=col.al[1],lwd=1.7)
legend('topright',bty='n',pch=c(16,16),col=col.al[c(1,2)],legend=c('Experiment','ADO'))
mtext('Young',side=3,cex=1.3,line=0.5)
polygon(c(seq(1,192),rev(seq(1,192))),
        c(klq[3,,1,2],rev(klq[1,,1,2])),
        col=col.a[2],border=F)
lines(seq(1,192),klq[2,,1,2],col=col.al[2],lwd=1.7)

plot(seq(1,192),ylim=c(0,max(klq[,,2,])),type='n',axes=F,ann=F)
box(bty='l')
axis(1,at=c(1,50,100,150,192))
axis(2,las=2)
mtext('Old',side=3,cex=1.3,line=0.5)
legend('topright',bty='n',pch=c(16,16),col=col.al[c(1,3)],legend=c('Experiment','ADO'))
polygon(c(seq(1,192),rev(seq(1,192))),
        c(klq[3,,2,1],rev(klq[1,,2,1])),
        col=col.a[1],border=F)
lines(seq(1,192),klq[2,,2,1],col=col.al[1],lwd=1.7)
polygon(c(seq(1,192),rev(seq(1,192))),
        c(klq[3,,2,2],rev(klq[1,,2,2])),
        col=col.a[3],border=F)
lines(seq(1,192),klq[2,,2,2],col=col.al[3],lwd=1.7)
mtext('Trial',side=1,outer=T,cex=1.3,line=1.2)
dev.off()

#### Figure 7: d' by participant non-contaminant model, (example used is elderly 14) ####
load("data/signaldt.RData")
ylable <- c('Criterion','Lure 1','Lure 2','Lure 3','Lure 4','Lure 5','New')
age.name <- c('Young','Elderly')
for(aa in 1:ages){
  for(pp in 1:n.sub){
    file.name <- paste(c('individualD_',bquote(.(age.name[aa])),'_nc',bquote(.(pp)),'.pdf'),collapse='')
    pdf(file=paste(c('figures/individualdprime_noncontaminant/',file.name),collapse=''))
    par(oma=c(2,2,2,.7))
    layout(rbind(c(1,2),
                 c(3,4),
                 c(5,6),
                 c(7,8)))
    par(mai=c(0.3,0.5,0.2,0.1),xaxs="i",yaxs='i')
    counts <- c()
    for(ss in 1:7){
      counts <- append(counts,which(sdt$optimal$storder[,1,pp,aa]==ss)[-c(1:table(results$st[,pp,aa])[ss])])
      plot(0,0,xlim=c(-1,trials),type="n",axes=F,ann=F,
           ylim=c(min(sdt$data$ci[,1,,pp,aa],sdt$optimal$ci[,1,,pp,aa]),
                  max(sdt$data$ci[,2,,pp,aa],sdt$optimal$ci[,2,,pp,aa])))
      if(ss==7){
        mtext(expression(paste("d'")),side=2,las=2,line=2.3,cex=1.3)
        mtext("Trial",side=1,line=3,cex=1.5)
      }
      mtext(paste(ylable[ss]),side=3)
      axis(1,at=c(1,50,100,150,192),cex.axis=1.25)
      axis(2,las=2,cex.axis=1.25)
      polygon(x=c(seq(1,times),rev(seq(1,times))),border=F,col=col.a[1],
              y=c(sdt$data$ci[ss,2,,pp,aa],rev(sdt$data$ci[ss,1,,pp,aa])))
      lines(seq(1,times),sdt$optimal$ci[ss,1,,pp,aa],col=col.al[aa+1],lwd=1.2)
      lines(seq(1,times),sdt$optimal$ci[ss,2,,pp,aa],col=col.al[aa+1],lwd=1.2)
    }
    plot(0,0,ylim=c(min(sdt$data$ut[,pp,aa],sdt$optimal$ut[,pp,aa]),
                    max(sdt$data$ut[,pp,aa],sdt$optimal$ut[,pp,aa])),
         type="n",axes=F,ann=F,xlim=c(0,dim(sdt$data$ut)[1]))
    legend('topright',bty='n',col=c(col.al[1],col.al[aa+1]),
           legend=c('Experiment','ADO'),cex=1.3,
           pch=c(16,16))
    lines(smooth.spline(seq(1,dim(sdt$data$ut)[1]),sdt$data$ut[,pp,aa]),col=col.al[1],lwd=2)
    lines(smooth.spline(seq(1,dim(sdt$optimal$ut)[1]),sdt$optimal$ut[,pp,aa]),col=col.al[aa+1],lwd=2)
    points(seq(1,dim(sdt$data$ut)[1]),sdt$data$ut[,pp,aa],col=col.a[1],pch=16,cex=0.8)
    points(seq(1,dim(sdt$optimal$ut)[1]),sdt$optimal$ut[,pp,aa],col=col.a[aa+1],pch=16,cex=0.8)
    axis(1,at=c(1,50,100,150,192),cex.axis=1.25)
    axis(2,las=2,cex.axis=1.25,hadj = 0.8)
    mtext("KL Divergence",side=2,line=2.7,cex=1.1)
    abline(v=min(counts),col='slategray',lty=2)
    dev.off()
  }
}
#### Figure 8: Stimulus used by method, examples used is young 15 (fig. 7) ####
load('data/contaminant_sdt.RData')
ylable <- c('Old','Lure 1','Lure 2','Lure 3','Lure 4','Lure 5','New')
age.name <- c('Young','Elderly')
for(aa in 1:ages){
  for(pp in 1:n.sub){
    file.name <- paste(c('storder',bquote(.(age.name[aa])),'_',bquote(.(pp)),'.pdf'),collapse='')
    pdf(file=paste(c('figures/indstimulusorder/',file.name),collapse=''))
    par(oma=c(2,2,2,.7))
    layout(rbind(c(1,2),
                 c(3,4),
                 c(5,6),
                 c(7,8)))
    par(mai=c(0.3,0.5,0.2,0.1),xaxs="i",yaxs='i')
    counts <- c()
    for(ss in 1:7){
      counts <- append(counts,which(sdt$optimal$storder[,1,pp,aa]==ss)[-c(1:table(results$st[,pp,aa])[ss])])
      plot(0,0,xlim=c(-1,trials),type="n",axes=F,ann=F,
           ylim=c(-0.05,1.2))
      if(ss==7){
        mtext("Trial",side=1,line=3,cex=1.5)
      }
      mtext(paste(ylable[ss]),side=3)
      box(bty='l')
      axis(1,at=c(1,50,100,150,192),cex.axis=1.25)
      axis(2,at=seq(0,1,0.2),labels = c('0',seq(0.2,0.8,0.2),'1'),las=2,cex.axis=1.25)
      optimal.order <- results$st[sdt$optimal$expdes[,pp,aa],pp,aa]
      lines(seq(1,times),cumsum(optimal.order==ss)/sum(optimal.order==ss),
            col=col.al[aa+1],lwd=1.2)
      points(which(sdt$optimal$storder[,1,pp,aa]==ss),
             c(cumsum(optimal.order==ss)/sum(optimal.order==ss))[which(sdt$optimal$storder[,1,pp,aa]==ss)],
             pch=16,col=col.al[2])
      lines(seq(1,times),cumsum(results$st[,pp,aa]==ss)/sum(results$st[,pp,aa]==ss),
            col=col.al[1],lwd=1.2)
    }
    mtext('Proportion of Stimulus',side=2,outer=T,line=-0.8,cex=1.5)
    dev.off()
  }
}

#### Figure 9: proportion of used stimulus by trial ####
load('data/contaminant_sdt.RData')
col.l <- c(rgb(0.5843,0.3216,0.3176),
           rgb(0.9176,0.8510,0.5451),
           rgb(0.8322,0.8071,0.5859),
           rgb(0.7467,0.7631,0.6267),
           rgb(0.6612,0.7192,0.6675),
           rgb(0.5757,0.6753,0.7082),
           rgb(0.4902,0.6314,0.7490))
prop.time.p <- array(NA,dim=c(7,192,20,2))

for(a in 1:ages){
  for(p in 1:n.sub){
    pp <- table(results$st[,p,a])
    for(t in 1:trials){
      for(s in 1:parameters){
        prop.time.p[s,t,p,a] <- sum(results$st[sdt$optimal$expdes[1:t,p,a],p,a]==s)/pp[s]
      }
    } 
  }
}

age.name <- c('Young','Elderly')
ylable <- c('Old','Lure 1','Lure 2','Lure 3','Lure 4','Lure 5','New')
dist <- c(rep(c(1,2),3),1)
pdf(file='figures/proportionst_time.pdf',width=8,height = 4)
par(oma=c(2.3,3,1,1))
layout(t(c(1,2)))
par(mai=c(0.3,0.3,0.2,0.1),xaxs="i",yaxs='i')
for(a in 1:ages){
  plot(0,0,xlim=c(0,trials+1),ylim=c(0,1.05),axes=F,ann=F,type='n')
  box(bty='l')
  axis(1,at=c(1,50,100,150,192),cex.axis=1.2)
  mtext(paste(age.name[a]),side=3,cex=1.3,line=0.5)
  if(a==1){
    axis(2,at=seq(0,1,0.2),labels=c('0',seq(0.2,0.8,0.2),'1'),cex.axis=1.2,las=2,
         hadj = 0.8)
  }
  else{
    axis(2,at=seq(0,1,0.2),labels=rep('',6),cex.axis=1.2,las=2,
         hadj = 0.8)
    legend('bottomright',legend = ylable,col=col.l,lty=dist,lwd=2,bty='n',cex=0.86)
  }
  for(s in 1:parameters){
    lines(seq(1,trials),rowMeans(prop.time.p[s,,,a]),col=col.l[s],type='s',
          lty=dist[s],lwd=3)
  }
}
mtext(text = 'Proportion of stimulus presented',cex=1.5,outer=T,side=2,line=1.3)
mtext('Trial',side=1,outer=T,cex=1.5,line=1.2)
dev.off()

#### Figure 10: KL for simulated experimental designs and ado ####
load('data/simulation_multiple_designs.Rdata')
pdf(file=paste(c('figures/','simulations_KL.pdf'),collapse=''),width = 8)
par(oma=c(3,4,2,.7))
layout(matrix(seq(1,6),nrow = 2,ncol = 3,byrow = T))
par(mai=c(0.2,0.2,0.2,0.1),xaxs="i",yaxs='i')
for(aa in 1:2){
  for(pp in 1:3){
    plot(0,0,ylim=c(min(sdt$data$ut[,,,],sdt$optimal$ut[,pp,aa]),
                    max(sdt$data$ut[,,,],sdt$optimal$ut[,pp,aa])),
         type="n",axes=F,ann=F,xlim=c(0,dim(sdt$data$ut)[1]))
    box(bty='l')
    for(dd in 1:5){
      lines(smooth.spline(seq(1,dim(sdt$data$ut)[1]),sdt$data$ut[,dd,pp,aa]),col=col.al[1],lwd=2)
      points(seq(1,dim(sdt$data$ut)[1]),sdt$data$ut[,dd,pp,aa],col=col.a[1],pch=16,cex=0.8)
    }
    if(pp==3){
      legend('topright',bty='n',col=c(col.al[1],col.al[aa+1]),
             legend=c('Experiment','ADO'),cex=1.8,
             pch=c(16,16))
    }
    lines(smooth.spline(seq(1,dim(sdt$optimal$ut)[1]),sdt$optimal$ut[,pp,aa]),col=col.al[aa+1],lwd=2)
    points(seq(1,dim(sdt$optimal$ut)[1]),sdt$optimal$ut[,pp,aa],col=col.a[aa+1],pch=16,cex=0.8)
    axis(1,at=c(1,50,100,150,192),cex.axis=1.5)
    axis(2,las=2,cex.axis=1.5,hadj = 0.85)
    if(pp==1&aa==2){
      mtext("KL Divergence",side=2,line=1.9,outer=T,cex=1.5)
      mtext("Trial",side=1,line=1.8,outer=T,cex=1.5)
    }
  }
}
dev.off()

#### Figure 11 was done in MATLAB, code is in 'src/draftFigure.m' ####
#### Appendix ####
# Plots 1 and 2 KL by participant and eage group
ylable <- c('Criterion','Lure 1','Lure 2','Lure 3','Lure 4','Lure 5','New')
age.name <- c('Young','Elderly')
for(aa in 1:ages){
    file.name <- paste(c('KL_',bquote(.(age.name[aa])),'.pdf'),collapse='')
    pdf(file=paste(c('figures/appendix/',file.name),collapse=''),width = 8)
    par(oma=c(3,3,2,.7))
    layout(matrix(seq(1,20),nrow = 4,ncol = 5,byrow = T))
    par(mai=c(0.2,0.2,0.2,0.1),xaxs="i",yaxs='i')
  for(pp in 1:n.sub){
    counts <- c()
    for(ss in 1:7){
      counts <- append(counts,which(sdt$optimal$storder[,1,pp,aa]==ss)[-c(1:table(results$st[,pp,aa])[ss])])
    }
    plot(0,0,ylim=c(min(sdt$data$ut[,pp,aa],sdt$optimal$ut[,pp,aa]),
                    max(sdt$data$ut[,pp,aa],sdt$optimal$ut[,pp,aa])),
         type="n",axes=F,ann=F,xlim=c(0,dim(sdt$data$ut)[1]))
    if(pp==20){
      legend('topright',bty='n',col=c(col.al[1],col.al[aa+1]),
           legend=c('Experiment','ADO'),cex=1.3,
           pch=c(16,16))
      }
    lines(smooth.spline(seq(1,dim(sdt$data$ut)[1]),sdt$data$ut[,pp,aa]),col=col.al[1],lwd=2)
    lines(smooth.spline(seq(1,dim(sdt$optimal$ut)[1]),sdt$optimal$ut[,pp,aa]),col=col.al[aa+1],lwd=2)
    points(seq(1,dim(sdt$data$ut)[1]),sdt$data$ut[,pp,aa],col=col.a[1],pch=16,cex=0.8)
    points(seq(1,dim(sdt$optimal$ut)[1]),sdt$optimal$ut[,pp,aa],col=col.a[aa+1],pch=16,cex=0.8)
    axis(1,at=c(1,100,192),cex.axis=1.4)
    axis(2,las=2,cex.axis=1.4)
    if(pp==11){
      mtext("KL Divergence",side=2,line=3,cex=1.5)
    }
    if(pp==18){
      mtext("Trial",side=1,line=3,cex=1.5)
    }
    abline(v=min(counts),col='slategray',lty=2)
  }
  dev.off()
}
#### Revisions ####
# Plot 1: number of trials to proportion of KL ----
load('data/contaminant_sdt.RData')
max_kl <- matrix(NA, nrow = n.sub, ncol = ages)
max_kl_young_optimal <- apply(X = sdt$optimal$ut[,,1], MARGIN = 2, FUN = max)
max_kl_young_data <- apply(X = sdt$data$ut[,,1], MARGIN = 2, FUN = max)

max_kl[,1] <- apply(X = cbind(max_kl_young_data, max_kl_young_optimal),
                    MARGIN = 1, FUN = max)

max_kl_elder_optimal <- apply(X = sdt$optimal$ut[,,2], MARGIN = 2, FUN = max)
max_kl_elder_data <- apply(X = sdt$data$ut[,,2], MARGIN = 2, FUN = max)

max_kl[,2] <- apply(X = cbind(max_kl_elder_data, max_kl_elder_optimal),
                    MARGIN = 1, FUN = max)

percentage_below <- c(0.99,rev(seq(0.2,0.99,0.05)))

trial_lower_ut_optimal <- array(NA, dim = c(length(percentage_below),
                                            length(max_kl[,1]),
                                            ages))

trial_lower_ut_data <- array(NA, dim = c(length(percentage_below),
                                         length(max_kl[,1]),
                                         ages))

for(pp in 1:ages){
  for(i in 1:n.sub){
    count <- 0
    for(j in percentage_below){
      count <- count + 1
      trial_lower_ut_optimal[count,i,pp] <- ifelse(test = is.finite(
        min(which(sdt$optimal$ut[,i,pp] <= max_kl[i,pp]*j))),
        yes = min(which(sdt$optimal$ut[,i,1] <= max_kl[i,pp]*j)),
        no = 192)
      
      trial_lower_ut_data[count,i,pp] <- ifelse(test = is.finite(
        min(which(sdt$data$ut[,i,pp] <= max_kl[i,pp]*j))),
        yes = min(which(sdt$data$ut[,i,pp] <= max_kl[i,pp]*j)),
        no = 192)
    }
  }  
}

klq <- array(NA,dim=c(3,length(percentage_below),2,2))
# obtain the quantiles of the KL divergence for each age group and each trial
klq[1:3,,1,1] <- apply(X = trial_lower_ut_data[,,1], MARGIN = 1, 
                       FUN = quantile, p = c(0.25,0.5,0.75))
klq[1:3,,2,1] <- apply(X = trial_lower_ut_data[,,2], MARGIN = 1, 
                       FUN = quantile, p = c(0.25,0.5,0.75))
klq[1:3,,1,2] <- apply(X = trial_lower_ut_optimal[,,1], MARGIN = 1, 
                       FUN = quantile, p = c(0.25,0.5,0.75))
klq[1:3,,2,2] <- apply(X = trial_lower_ut_optimal[,,2], MARGIN = 1, 
                       FUN = quantile, p = c(0.25,0.5,0.75))

pdf(file='figures/trials-prop-KL.pdf',width=8,height = 4)

par(oma=c(2.3,3,1,1))
layout(t(c(1,2)))
par(mai=c(0.3,0.3,0.2,0.1),xaxs="i",yaxs='i')

plot(percentage_below, klq[1,,1,1] , ylim = c(0, max(klq[,,1,])), 
     xlim = rev(range(percentage_below)), type='n', axes = F, ann = F)

box(bty = 'l')
axis(1)
axis(2, las = 2)

polygon(c(percentage_below, rev(percentage_below)),
        c(klq[3,,1,1], rev(klq[1,,1,1])),
        col = col.a[1], border = FALSE)
lines(percentage_below, klq[2,,1,1], col = col.al[1], lwd = 1.7)

polygon(c(percentage_below, rev(percentage_below)),
        c(klq[3,,1,2], rev(klq[1,,1,2])),
        col = col.a[2], border = FALSE)
lines(percentage_below, klq[2,,1,2], col = col.al[2], lwd = 1.7)

mtext('Trials', side = 2, line = 2.5, cex = 1.3)
legend('topleft', bty = 'n',pch = c(16, 16), col = col.al[c(1,2)],
       legend = c('Experiment', 'ADO'))
mtext('Young', side = 3, cex = 1.3, line = 0.5)

plot(percentage_below, klq[1,,2,1], ylim = c(0, max(klq[,,2,])),
     xlim = rev(range(percentage_below)), type = 'n', axes = FALSE, ann = FALSE)

box(bty = 'l')
axis(1)
axis(2, las = 2)

polygon(c(percentage_below, rev(percentage_below)),
        c(klq[3,,2,1], rev(klq[1,,2,1])),
        col = col.a[1], border = FALSE)
lines(percentage_below, klq[2,,2,1], col = col.al[1], lwd = 1.7)

polygon(c(percentage_below, rev(percentage_below)),
        c(klq[3,,2,2], rev(klq[1,,2,2])),
        col = col.a[3], border = FALSE)
lines(percentage_below, klq[2,,2,2],col = col.al[3], lwd = 1.7)

mtext('Old',side = 3,cex = 1.3,line = 0.5)
legend('topleft', bty='n', pch=c(16,16), col = col.al[c(1,3)],
       legend = c('Experiment','ADO'))

mtext("Proportion of max KL", side = 1, outer = TRUE , cex = 1.3, line = 1.2)

dev.off()

# Plot 2: Savage dickey density ratio for boundary = 0.01 ----
load("data/bf-contaminant-model.RData")

pdf(file='figures/bayes-factors.pdf',width=8,height = 4)
par(oma=c(2.3,4,1,1))
layout(t(c(1,2)))
par(mai=c(0.3,0.3,0.2,0.1))

plot(seq(1, 20), -log(bf[,1]), axes = FALSE, ann = FALSE, 
     xlim = c(1,20), pch = 21, col = rep(col.al[2], each = n.sub),
     bg = rep(col.a[2], each = n.sub),
     ylim = c(-4,4))
abline(h = 0, lty = 3, col = "#36454F", lwd = 1.2)

box(bty = "l")
axis(1, at = c(seq(1,20,4),20))
axis(2, las = 2, at = log(c(sqrt(10), 10, sqrt(10)^3)), 
     labels = c(expression(10^(1/2)),10,expression(10^(3/2))), hadj = 0.84)
axis(2, las = 2, at = c(-log(sqrt(10)), -log(10), -log(sqrt(10)^3), 0), 
     labels = c(expression(1/10^(1/2)), "1/10", expression(1/10^(3/2)), 1),
     hadj = 0.885)

mtext('Young', side = 3, cex = 1.3, line = 0.5)
mtext('Bayes Factor', side = 2, line = 4.3, cex = 1.3)

plot(seq(1, 20), -log(bf[,2]), axes = FALSE, ann = FALSE, 
     xlim = c(1,20), pch = 21, col = rep(col.al[3], each = n.sub),
     bg = rep(col.a[3], each = n.sub),
     ylim = c(-4,4))
abline(h = 0, lty = 3, col = "#36454F", lwd = 1.2)

box(bty = "l")
axis(1, at = c(seq(1,20,4),20))
axis(2, las = 2, at = log(c(sqrt(10), 10, sqrt(10)^3)), 
     labels = rep("", 3))
axis(2, las = 2, at = c(-log(sqrt(10)), -log(10), -log(sqrt(10)^3), 0), 
     labels = rep("", 4))

mtext('Old',side = 3,cex = 1.3,line = 0.5)

mtext("Participant", side = 1, outer = TRUE , cex = 1.3, line = 1.2)

dev.off()
