# function to create a scatter plot of two series,
# with circles as symbols and simple lines as error-bars
# includes horizontal lines for averages and
# complext y-label

# data is stored in S at the end of this script

mk_fig5<- function(out='screen',S=data_fig5){
  
  # N and C values are now uptake over 2 days,
  # to express in values per day they are divided by 2
  S$C <- S$C/2
  S$N <- S$N/2
  
  m_C <- aggregate(S$C~S$G,FUN=mean)
  s_C <- aggregate(S$C~S$G,FUN=sd)
  se_C <- aggregate(S$C~S$G,FUN=function(x) sd(x)/sqrt(length(x)) )
  mC <- data.frame(g=m_C[,1], mean=m_C[,2], se=se_C[,2]) 
  
  m_N <- aggregate(S$N~S$G,FUN=mean)
  s_N <- aggregate(S$N~S$G,FUN=sd)
  se_N <- aggregate(S$N~S$G,FUN=function(x) sd(x)/sqrt(length(x)) )
  mN <- data.frame(g=m_N[,1], mean=m_N[,2], se=se_N[,2]) 
  
  # define vectors for location on x-axis and vertical 
  # (standard error) lines
  
  xf <- c(1:5)
  xn <- c(1:5)
  xfl <- c( rbind(xf,xf,rep(NA,5)) )
  yfl <- c( rbind( mC$mean+mC$se, mC$mean-mC$se, rep(NA,5)) )
  xnfl <- c( rbind(xn,xn,rep(NA,5)) )
  ynfl <- c( rbind( mN$mean+mN$se, mN$mean-mN$se, rep(NA,5)) )
  
  Cavg <- mean(mC$mean)
  Navg <- mean(mN$mean)
  
  label_yaxis1 <- 'Amino acid assimilation'
  label_yaxis2 <- expression('('*μmol~C~or~N[~tracer]~mmol~C~or~N[~sponge]^{-1}~d^{−1}*')')
  label_yaxis3 <- expression('('*μmol~C~or~N[~tracer]~(mmol~C~or~N[~sponge]~d)^{−1}*')')
  
  
  if(out=='png'){ png('fig5.png') }
  if(out=='pdf'){ pdf('fig5.pdf') }
  
  # svg doesn't work with expression
  # if(out=='svg'){ svg('fig5.svg') }

  op <- par(no.readonly = TRUE)
  par(mar=c(4,7,0.5,0.5))
  
  plot(xf, mC$mean, type='n', xaxt="n", bty='n', las=1,
       xlim = c(0, 6), ylim = c(0, 3), xlab = "", xaxs="i",
       ylab = '', cex.lab=1.2, cex.axis=1.4)
  
  lines(c(-0.1,5.2),c(Cavg,Cavg),col='grey',lty=2,lwd=2)
  lines(c(-0.1,5.2),c(Navg,Navg),col='grey',lty=2,lwd=2)
  lines(xfl, yfl, col = "black")
  lines(xnfl, ynfl, col = "black")
  
  points(xf, mC$mean, col = "black", pch = 21, bg = "grey", cex = 2)
  points(xn, mN$mean, col = "black", pch = 21, bg = "black", cex = 2)
  
  # mtext( paste(S$G, 'g'), at=1:5, side=1, line=0.3, cex=1.2)
  # mtext( 'uptake (umol/day)', at = 0, side=3, line=0, cex=1.3)
  # axis(1, at = 0:5, labels = c('', paste( unique(S$G),'g')), cex.lab=1.4, cex.axis=1.4, pos=0)
  axis(1, at = 0:5, labels = c('', unique(S$G)), cex.lab=1.4, cex.axis=1.4, pos=0)
  title(ylab = label_yaxis1, cex.lab = 1.4, line = 5.7)
  title(ylab = label_yaxis2, cex.lab = 1.4, line = 3.3)
  # alternative 2nd line in y-axis label
  # title(ylab = label_yaxis3, cex.lab = 1.4, line = 3.3)
  
  title(xlab = expression(italic(g)*'-force'), cex.lab = 1.4, line = 1.8)
  
  
  text(5.5, Cavg, "C", cex = 1.5, font = 1, adj = 0, col="darkgrey")
  text(5.5, Navg, "N", cex = 1.5, font = 1, adj = 0)
  
  par( op )
  
  if(out!='screen'){ dev.off() }
  
}



# data to be plotted

data_fig5 <- structure(list( osc = c(1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L), 
                     G = c(1.0, 2.5, 5.0, 10.0, 20.0, 1.0, 2.5, 5.0, 10.0, 20.0), 
                     C = c(0.997897694, 0.648214747, 0.974571948, 0.662205429, 0.711240641, 0.837062031, 0.782400568, 0.87364085, 0.845282945, 0.725389001), 
                     N = c(4.933197895, 4.127039988, 5.851914928, 3.291699288, 2.912717366, 4.483257766, 5.013676084, 5.467241384, 4.420443553, 2.972087615)), 
               .Names = c("osc","G", "C", "N"), class = "data.frame", row.names = c(NA, -10L))
