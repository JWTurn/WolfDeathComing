library("ggplot2")


####looking at rss of selecting an area x distance from road, when 100 m from roads, up to 10 km from roads
###(1+Δhi/hi)^Bi
###keep the change in distance the same but change the start distance


hj = c(250:750)
dhj = 250
ehj=exp(hj)

####Model 1 selecting roads

bjday = 0.3287 
bjday_lo = 0.2917
bjday_hi = 0.3713

bjtwi = 0.3779
bjtwi_lo = 0.3327
bjtwi_hi = 0.4298
  
bjnight = 0.2214  	
bjnight_lo = 0.1922 
bjnight_hi = 0.2529
  
rssroadd = (1+(hj/(hj-dhj)))^bjday
rssroadd_lo = (1+(hj/(hj-dhj)))^bjday_lo
rssroadd_hi = (1+(hj/(hj-dhj)))^bjday_hi

rssroadt = (1+(hj/(hj-dhj)))^bjtwi
rssroadt_lo = (1+(hj/(hj-dhj)))^bjtwi_lo
rssroadt_hi = (1+(hj/(hj-dhj)))^bjtwi_hi

rssroadn = (1+(hj/(hj-dhj)))^bjnight
rssroadn_lo = (1+(hj/(hj-dhj)))^bjnight_lo
rssroadn_hi = (1+(hj/(hj-dhj)))^bjnight_hi

r =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7)
r = r + geom_line(aes(x=(hj-250),y=log(rssroadd), colour = "Day"), size = 1) 
r = r + geom_line(aes(x=(hj-250),y=log(rssroadd_lo), colour = "Day"), size = 1, lty = 3) 
r = r + geom_line(aes(x=(hj-250),y=log(rssroadd_hi), colour = "Day"), size = 1, lty = 3) 
r = r + geom_line(aes(x=(hj-250),y=log(rssroadt), colour = "Twilight"), size = 1) 
r = r + geom_line(aes(x=(hj-250),y=log(rssroadt_lo), colour = "Twilight"), size = 1, lty = 3) 
r = r + geom_line(aes(x=(hj-250),y=log(rssroadt_hi), colour = "Twilight"), size = 1, lty = 3) 
r = r + geom_line(aes(x=(hj-250),y=log(rssroadn), colour = "Night"),size = 1) 
r = r + geom_line(aes(x=(hj-250),y=log(rssroadn_lo), colour = "Night"),size = 1, lty = 3) 
r = r + geom_line(aes(x=(hj-250),y=log(rssroadn_hi), colour = "Night"),size = 1, lty = 3) 
r = r + theme_bw()  + theme(
  #panel.background =element_rect(colour = "black", fill=NA, size=1),
  panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black", size = .7))
r = r + theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20))
r = r + theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
              axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) 
r = r + ylab("ln(RSS)") + xlab("Distance to Road (m)")
r = r + ylim(-0.01,1.5)
r = r +scale_colour_manual("", 
                             values = c("gray", "black", "gray33"))  
r = r +  theme(legend.key = element_blank()) + theme(legend.position = c(.3,.9)) + theme(legend.text = element_text(size = 20))

print(r)




####road effect distance 
xd = 500/(1-(2^(-1/bjday)))

xd - 250

xn = 500/(1-(2^(-1/bjnight)))

xn - 250

xt = 500/(1-(2^(-1/bjtwi)))

xt - 250


### rss of selecting a value of NDVI given you are x from road
###exp[Δhi(βi + 2iβi2 + Δhiβhi2 + hjβij)]
###make Δhi on number and change hj

####for day with changing road distance
###road distance is y
###NDVI is x, normalized values range from 0 to 100, delta x is dx

hix1= 50
hix2 = 20

############ day is 0.1128
bid = 0.1128
bid_lo = 0.0863
bid_hi = 0.1407

### quadratic for ndvi is -0.0006 ###covarcalc
bi2d = -0.000647
bi2d_lo = -0.000524
bi2d_hi = -0.000776

###coef for road:NDVI is ###covarcalc
bij =  -0.008954
bij_lo = -0.006335
bij_hi = -0.011712


hj =seq(0,750)
lnhj = log(hj+1)

###equation

rssndvi = (hix1 - hix2) * ((bid+(bi2d*(hix1-hix2))) + (bij*(lnhj)))

rssndvi_lo = (hix1 - hix2) * ((bid_lo+(bi2d_lo*(hix1-hix2))) + (bij_lo*(lnhj)))
rssndvi_hi = (hix1 - hix2) * ((bid_hi+(bi2d_hi*(hix1-hix2))) + (bij_hi*(lnhj)))

n2 = ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7)
n2 = n2 + geom_line(aes(x=hj,y=log(rssndvi)), colour = "black", size = 1,lty=1) 
n2 = n2+ geom_line(aes(x=hj,y=log(rssndvi_lo)), colour = "black", size = 1,lty=3) 
n2 = n2+ geom_line(aes(x=hj,y=log(rssndvi_hi)), colour = "black", size = 1,lty=3) 
n2 = n2+  theme_bw()  + theme(
  #panel.background =element_rect(colour = "black", fill=NA, size=1),
  panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black", size = .7))
n2 = n2 +theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20))
n2 = n2 + theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
                axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) 

n2 = n2 + xlab(expression('Distance to Road (m) at x'[1])) + ylab(" ") + ggtitle("(d)")  
print(n2)
 



multiplot(n,xr, rd,n2, cols=2)

xr


