##
##  g a m m a z . R  Test suite
##


gammaz <- pracma::gammaz

##  Problematic on Solaris (2012-01-25)
# y <- seq(from=0,to=5,by=0.5)
# # z0 <- lngamma_complex(1+y*1i)
# z0 <- c(0.0000000+0.0000000i, -0.1909455-0.2440583i, -0.6509232-0.3016403i,
#        -1.2344831-0.1629398i, -1.8760788+0.1296463i, -2.5499068+0.5426044i,
#        -3.2441443+1.0533508i, -3.9524671+1.6461926i, -4.6710996+2.3096981i,
#        -5.3976062+3.0351970i, -6.1303241-2.4672867i)
# 
# all.equal(gammaz(1+y*1i), exp(z0), tolerance = 1e-7)
