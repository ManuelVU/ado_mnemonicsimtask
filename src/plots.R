# This code can be used to generate the figures in the "Adaptive Design Optimization for a Mnemonic Similarity Task"
# Figures in the paper are numbered in the code, however, it does not generate particular examples but figures for 
# all participants when needed. Figures 2 and 3 are examples and are generated independently.

#### Load data and libraries needed ####
library(R.matlab)
load("data/memory.RData")
load('data/contaminant_sdt.RData')
parameters <- 7
times <- 192
trials <- 192
n.sub <- 20
ages <- 2

# load color pallete
colors_pantone <- readMat('figures/PantoneFall2016.mat')
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
        mtext(expression(paste("d'")),side=2,las=2,line=2.3)
        mtext("Trials",side=1,line=2.4)
      }
      mtext(paste(ylable[ss]),side=3)
      axis(1,at=c(1,50,100,150,192),cex.axis=1.15)
      axis(2,las=2,cex.axis=1.15)
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
    axis(1,at=c(1,50,100,150,192),cex.axis=1.15)
    axis(2,las=2,cex.axis=1.15)
    mtext("KL Divergence",side=2,line=2.5)
    abline(v=min(counts),col='slategray',lty=2)
    dev.off()
  }
}

#### Figure 6: KL quantiles by age group ####
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

#### Figure 7: Stimulus used by method, examples used aisre young 15 (fig. 7) ####
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
        mtext("Trials",side=1,line=3,cex=1.5)
      }
      mtext(paste(ylable[ss]),side=3)
      box(bty='l')
      axis(1,at=c(1,50,100,150,192),cex.axis=1.15)
      axis(2,at=seq(0,1,0.2),labels = c('0',seq(0.2,0.8,0.2),'1'),las=2,cex.axis=1.15)
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

#### FIgure 8: KL for simulated experimental designs and ado ####
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
      mtext("Trials",side=1,line=1.8,outer=T,cex=1.5)
    }
  }
}
dev.off()

#### Figure 9 was done in matlab, code is in 'src/draftFigure.m' ####

#### Appendinx ####
# Plots 1 and 2 KL by participant and eage group
ylable <- c('Criterion','Lure 1','Lure 2','Lure 3','Lure 4','Lure 5','New')
age.name <- c('Young','Old')
for(aa in 1:ages){
    file.name <- paste(c('KL_',bquote(.(age.name[aa])),'.pdf'),collapse='')
    pdf(file=paste(c('Dropbox/signal_detection/smp2021/figures/KLage/',file.name),collapse=''),width = 8)
    par(oma=c(2,2.3,2,.7))
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
    axis(1,at=c(1,50,100,150,192),cex.axis=1.15)
    axis(2,las=2,cex.axis=1.15)
    if(pp==11){
      mtext("KL Divergence",side=2,line=2.5)
    }
    if(pp==18){
      mtext("Trials",side=1,line=2.4)
    }
    abline(v=min(counts),col='slategray',lty=2)
  }
  dev.off()
}

# Plot 3 individual d' by participant (example used is participant young 15)
load('data/signaldt.RData')
ylable <- c('Criterion','Lure 1','Lure 2','Lure 3','Lure 4','Lure 5','New')
age.name <- c('Young','Elderly')
for(aa in 1:ages){
  for(pp in 1:n.sub){
    file.name <- paste(c('individualD_',bquote(.(age.name[aa])),'_nc',bquote(.(pp)),'.pdf'),collapse='')
    pdf(file=paste(c('figures/appendix/individualdprime_noncontaminant/',file.name),collapse=''))
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
        mtext(expression(paste("d'")),side=2,las=2,line=2.3)
        mtext("Trials",side=1,line=2.4)
      }
      mtext(paste(ylable[ss]),side=3)
      axis(1,at=c(1,50,100,150,192),cex.axis=1.15)
      axis(2,las=2,cex.axis=1.15)
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
    axis(1,at=c(1,50,100,150,192),cex.axis=1.15)
    axis(2,las=2,cex.axis=1.15)
    mtext("KL Divergence",side=2,line=2.5)
    abline(v=min(counts),col='slategray',lty=2)
    dev.off()
  }
}
