# Plot giuh histogram and compare with TauDEM-based (Python) giuh
# TauDEM-based GIUH (Python) are provided for catchment IDs 2943, 3789, 4366, and 5232

ncat <- 1 # 1, 50, 80, 120
print (giuh_compute['divide_id'], n=ncat)

# Visualization
giuh_dist <- fromJSON(giuh_compute["fun.giuh_minute"]$fun.giuh_minute[ncat])
giuh_dist

giuh_R_ordinates = giuh_dist$frequency[1:6]  # ordinates
giuh_R_times = giuh_dist$v[1:6]              # times

# TWI-Py (Luciana's script) 
# cat 2943
giuh_py_times = seq(from=60, to=360,by=60)

if (ncat == 1) {
  giuh_py_ordinates = c(0.58,0.29,0.13,0.0,0.0,0.0)  
  cat = "2943"
}
  

# cat - 3789 (ncat=50)
if (ncat == 50) {
  giuh_py_ordinates = c(0.41,0.37,0.22,0.0,0,0)
  cat = "3789"
}
  

# cat-4366 (n=80)
if (ncat == 80) {
  giuh_py_ordinates = c(0.62,0.28,0.10,0.0,0.0,0.0)
  cat = "4366"
}
  

# cat-5232 (ncat=120)
if (ncat == 120) {
  giuh_py_ordinates = c(0.66,0.34,0.0,0.0,0.0,0.0)
  cat = "5232"
}
  
  
# to plot multiple bar plot in one figure
barplot(rbind(giuh_py_ordinates, giuh_R),beside=TRUE, col= c('red','black'))

df1 <- data.frame(giuh_py_times, giuh_R_ordinates, giuh_py_ordinates)
df2 <- reshape::melt(df1, id=c('giuh_py_times'))
dat_giuh <- c(giuh_R_times, giuh_py_times)

options(repr.plot.width = 16, repr.plot.height =10) 
ggplot(df2, aes(x=dat_giuh, y=value, fill=variable)) +
  geom_bar(stat='identity', position='dodge', width=10.5) +
  scale_fill_manual(name = "GIUH", 
                    label = c("R","Py"), 
                    values = c("red","black"))+
  theme(aspect.ratio = 0.5/1, plot.title = element_text(hjust = 0.5)) +
  labs(title = glue("GIUH comparison (cat-{cat})") , x = "Time [minutes]", y = "GIUH ordinates") +
  scale_x_continuous(limits = c(0, 360), breaks = c(0,60,120,180,240,300,360))
