
library("ggplot2")


###looking at the  odds of selecting NDVI 

###NDVI is x, normalized values range from 0 to 100, delta x is dx
hix1=seq(00,100, by=.4)
hix2 = 20


#####NDVI for times of day ssf 2 ###from narrow boot weighted

############ day is 0.1128
bid = 0.1128
bid_lo = 0.0863
bid_hi = 0.1407
  
### quadratic for ndvi is -0.0006 ###covarcalc
bi2d = -0.000647
bi2d_lo = -0.000524
bi2d_hi = -0.000776

####for day with changing road distance
###road distance is y
###coef for road:NDVI is ###covarcalc
bij =  -0.008954
bij_lo = -0.006335
bij_hi = -0.011712

####log(1), log(250+1), log(1000+1)

hj0 = log(1+1)
hj1 = log(10+1)
hj2= log(100+1)
hj3 = log(1000+1)


###odds for NDVI eqn with interaction with roads, using daytime
##exp[(Bi+Bi2(2hi+Δhi)+Bijhj)Δhi)

rssndvi0 = (hix1 - hix2) * ((bid+(bi2d*(hix1-hix2))) + (bij*(hj0)))
rssndvi0_lo = (hix1 - hix2) *((bid_lo+(bi2d_lo*(hix1-hix2)))+ (bij_lo*(hj0)))
rssndvi0_hi = (hix1 - hix2) *((bid_hi+(bi2d_hi*(hix1-hix2)))+ (bij_hi*(hj0)))


rssndvi1 = (hix1 - hix2) * ((bid+(bi2d*(hix1-hix2))) + (bij*(hj1)))
rssndvi1_lo = (hix1 - hix2) *((bid_lo+(bi2d_lo*(hix1-hix2)))+ (bij_lo*(hj1)))
rssndvi1_hi = (hix1 - hix2) *((bid_hi+(bi2d_hi*(hix1-hix2)))+ (bij_hi*(hj1)))

rssndvi2 = (hix1 - hix2) * ((bid+(bi2d*(hix1-hix2))) + (bij*(hj2)))
rssndvi2_lo = (hix1 - hix2) *((bid_lo+(bi2d_lo*(hix1-hix2)))+ (bij_lo*(hj2)))
rssndvi2_hi = (hix1 - hix2) *((bid_hi+(bi2d_hi*(hix1-hix2)))+ (bij_hi*(hj2)))

rssndvi3 = (hix1 - hix2) * ((bid+(bi2d*(hix1-hix2))) + (bij*(hj3)))
rssndvi3_lo = (hix1 - hix2) *((bid_lo+(bi2d_lo*(hix1-hix2)))+ (bij_lo*(hj3)))
rssndvi3_hi = (hix1 - hix2) *((bid_hi+(bi2d_hi*(hix1-hix2)))+ (bij_hi*(hj3)))

rd = ggplot() + geom_line(aes(x=hix1,y=rssndvi0, colour = "1 m"), size = 1)
rd = rd + geom_line(aes(x=hix1,y=rssndvi0_lo, colour = "1 m"), size = 1, lty = 3) 
rd = rd + geom_line(aes(x=hix1,y=rssndvi0_hi, colour = "1 m"), size = 1, lty = 3) 
#rd = rd + geom_line(aes(x=hix1,y=rssndvi1,colour = "10 m"), size = 1) 
#rd = rd + geom_line(aes(x=hix1,y=rssndvi1_lo, colour = "10 m"), size = 1, lty = 3) 
#rd = rd + geom_line(aes(x=hix1,y=rssndvi1_hi, colour = "10 m"), size = 1, lty = 3) 
#rd = rd + geom_line(aes(x=hix1,y=rssndvi2, colour = "100 m"),  size = 1) 
#rd = rd + geom_line(aes(x=hix1,y=rssndvi2_lo, colour = "100 m"), size = 1, lty = 3) 
#rd = rd + geom_line(aes(x=hix1,y=rssndvi2_hi, colour = "100 m"), size = 1, lty = 3)
rd = rd + geom_line(aes(x=hix1,y=rssndvi3, colour = "1000 m"),  size = 1) 
rd = rd + geom_line(aes(x=hix1,y=rssndvi3_lo, colour = "1000 m"), size = 1, lty = 3) 
rd = rd + geom_line(aes(x=hix1,y=rssndvi3_hi, colour = "1000 m"), size = 1, lty = 3)
rd = rd + geom_hline(yintercept = 0,colour = "black", size = 1, lty = 2)
#rd = rd + ylim(-1,1) 
rd = rd + theme_bw()  + theme(
  #panel.background =element_rect(colour = "black", fill=NA, size=1),
  panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black", size = .7))
rd = rd +theme(plot.title=element_text(size=20,hjust = 0.05),plot.title=element_text(size=30), axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20))
rd = rd + theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
                  axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) 
rd = rd  + xlab(expression('NDVI at x'[1])) + ylab(" ") + ggtitle("(b)") 

rd = rd +scale_colour_manual("", 
                             values = c("black", "gray"))  
rd = rd +  theme(legend.key = element_blank()) + theme(legend.position = c(.2,.95)) + theme(legend.text = element_text(size = 20))
print(rd)



