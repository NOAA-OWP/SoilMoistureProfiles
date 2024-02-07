# the scripts simply barplot and compare Whitebox- and TauDEM-based TWIs

library(ggpubr)
# TWI-Py (TauDEM-based Luciana's script)
py_freq = c(0.206146, 0.0, 0.0, 0.0, 0.0, 0.00128, 0.002561, 0.0, 0.002561, 0.002561, 0.005122, 0.003841, 0.002561, 0.0, 0.002561, 0.002561, 0.005122, 0.011524, 0.015365, 0.025608, 0.035851, 0.025608, 0.056338, 0.075544, 0.103713, 0.111396, 0.12548, 0.103713, 0.056338, 0.016645)
py_v = c(18.578403, 18.121514, 17.664624, 17.207735, 16.750845, 16.293955, 15.837066, 15.380176, 14.923287, 14.466397, 14.009507, 13.552618, 13.095728, 12.638838, 12.181949, 11.725059, 11.26817, 10.81128, 10.35439, 9.897501, 9.440611, 8.983722, 8.526832, 8.069942, 7.613053, 7.156163, 6.699274, 6.242384, 5.785494, 5.328605)

twi_R <- fromJSON(twi["fun.twi"]$fun.twi[1])

barplot(rbind(py_v, rev(twi_R$v)),beside=TRUE, col= c('red','black'))

py_x = seq(from=1, to=30,by=1)

df1 <- data.frame(py_x, rev(twi_R$v), py_v)

df2 <- reshape::melt(df1, id=c('py_x'))
dat_twi <- c( rev(twi_R$frequency), py_freq)

options(repr.plot.width = 8, repr.plot.height =10) 
ggplot(df2, aes(x=dat_twi, y=value, fill=variable)) +
  geom_bar(stat='identity', position='dodge', width=0.005) +
  scale_fill_manual(name = "TWI", 
                    label = c("R","Py"), 
                    values = c("red","black"))+
  theme(aspect.ratio = 1/1) +
  labs(title = "TWI comparison", x = "area fraction", y = "TWI")


df11 <- data.frame(py_x, py_v, rev(twi_R$v), rev(twi_R$frequency))

options(repr.plot.width = 16, repr.plot.height =4) 
p1 <- ggplot(df11, aes(x=rev(twi_R$frequency), y=rev(twi_R$v))) +
  geom_bar(stat='identity', position='dodge', width=0.003,fill='black') +
  labs(title = "TWI-R", x = "area fraction", y = "TWI")+
  coord_cartesian(ylim=c(0,20),xlim=c(0,0.2))

p2 <- ggplot(df11, aes(x=py_freq, y=py_v)) +
  geom_bar(stat='identity', position='dodge', width=0.003,alpha=1.0, fill='black') +
  labs(title = "TWI-Py", x = "area fraction", y = "TWI")+
  coord_cartesian(ylim=c(0,20),xlim=c(0,0.2))

# later add pre-computed TWI
#p3 <- ggplot(df11, aes(x=twi_vrt$frequency, y=twi_vrt$v)) +
#  geom_bar(stat='identity', position='dodge', width=0.003,alpha=1.0,fill='black') +
#  labs(title = "TWI-HF", x = "area fraction", y = "TWI")+
#  coord_cartesian(ylim=c(0,20),xlim=c(0,0.2))

#ggarrange(p1, p2, p3,
#          labels = c("A", "B", "C"),
#          ncol = 3, nrow = 1)

ggarrange(p1, p2,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

