# Plot fractional area distribution
# The script barplot and compares TWI computed using Whitebox and TauDEM, also plots
# pre-computed TWI provided by the hydrofabric team as a parquet file

# TauDEM-based TWI (Python) are provided for catchment IDs 3789 and 3789
library(ggpubr)
ncat <- 1 # 1, 50 <-- 

twi_dist <- fromJSON(twi["fun.twi"]$fun.twi[ncat])

twi_R_value = rev(twi_dist$v)           # TWI values = ln(A/tanB)
twi_R_area = rev(twi_dist$frequency)    # distribution of area corresponding to ln(A/tanB)

#pre-computed TWI
twi_pre_comp_dist = fromJSON(twi_pre_computed["fun.twi"]$fun.twi[ncat])

twi_pre_comp_value = rev(twi_pre_comp_dist$v)           # TWI values = ln(A/tanB)
twi_pre_comp_area = rev(twi_pre_comp_dist$frequency)    # distribution of area corresponding to ln(A/tanB)


#barplot(rbind(twi_py_value, twi_R_value),beside=TRUE, col= c('red','black'))

################################################################################
# TauDEM-based TWI (TWI-Py)
# cat 2943 (ncat=1)
if (ncat == 1) {
  twi_py_area = c(0.206146, 0.0, 0.0, 0.0, 0.0, 0.00128, 0.002561, 0.0, 0.002561, 0.002561, 0.005122, 
                0.003841, 0.002561, 0.0, 0.002561, 0.002561, 0.005122, 0.011524, 0.015365, 0.025608, 
                0.035851, 0.025608, 0.056338, 0.075544, 0.103713, 0.111396, 0.12548, 0.103713, 0.056338, 0.016645)
  twi_py_value = c(18.578403, 18.121514, 17.664624, 17.207735, 16.750845, 16.293955, 15.837066, 15.380176, 
             14.923287, 14.466397, 14.009507, 13.552618, 13.095728, 12.638838, 12.181949, 11.725059, 
             11.26817, 10.81128, 10.35439, 9.897501, 9.440611, 8.983722, 8.526832, 8.069942, 7.613053, 
             7.156163, 6.699274, 6.242384, 5.785494, 5.328605)
  cat = "2943"
}

# cat - 3789 (ncat=50)
if (ncat == 50) {
  twi_py_area = c(0.069710,0.000415,0.000830,0.000830,0.002075,0.002075,0.004149,0.003734,0.005809,0.004149,
                0.006224,0.006639,0.007884,0.005809,0.009129,0.018672,0.021577,0.020332,0.025311,0.046473,
                0.072614,0.084232,0.107884,0.109129,0.109544,0.106224,0.078423,0.046473,0.017842,0.005809)
  twi_py_value = c(17.528700,17.132669,16.736638,16.340608,15.944577,15.548546,15.152516,14.756485,14.360454,
             13.964424,13.568393,13.172362,12.776332,12.380301,11.984270,11.588239,11.192209,10.796178,
             10.400147,10.004117,9.608086,9.212055,8.816025,8.419994,8.023963,7.627932,7.231902,6.835871,
             6.439840,6.043810)
  cat = "3789"
}

py_x = seq(from=1, to=30, by=1)
################################################################################

# PLOT 
plot_twi_vertically = FALSE

df1 <- data.frame(py_x, twi_R_value, twi_py_value, twi_pre_comp_value)
df2 <- reshape::melt(df1, id=c('py_x'))
dat_twi <- c(twi_R_area, twi_py_area, twi_pre_comp_area)

if (plot_twi_vertically) {
  options(repr.plot.width = 8, repr.plot.height =10) 
  ggplot(df2, aes(x=dat_twi, y=value, fill=variable)) +
    geom_bar(stat='identity', position='dodge', width=0.005, alpha=0.7) +
    scale_fill_manual(name = "Method", 
                      label = c("R","Py", "Pre-comp"), 
                      values = c("red","black","blue"))+
    theme(aspect.ratio = 1.0/1) +
    labs(title = "TWI comparison", x = "area fraction", y = "TWI")  
} else {
  options(repr.plot.width = 8, repr.plot.height =10) 
  ggplot(df2, aes(x=value, y=dat_twi, fill=variable)) +
    geom_bar(stat='identity', position='dodge', width=0.5, alpha=0.7) +
    scale_fill_manual(name = "Method", 
                      label = c("R","Py", "Pre-comp"), 
                      values = c("red","black","blue"))+
    theme(aspect.ratio = 1.0/1) +
    labs(title = "TWI comparison", x = "TWI", y = "Area fraction")
}


# area fraction vertically




# Three subplots
df11 <- data.frame(py_x, twi_py_value, twi_R_area, twi_R_value, twi_pre_comp_value)


options(repr.plot.width = 16, repr.plot.height =10) 
p1 <- ggplot(df11, aes(x=twi_R_area, y=twi_R_value)) +
  geom_bar(stat='identity', position='dodge', width=0.003,fill='black') +
  labs(title = "R", x = "area fraction", y = "TWI")+
  coord_cartesian(ylim=c(0,20),xlim=c(0,0.2))

p2 <- ggplot(df11, aes(x=twi_py_area, y=twi_py_value)) +
  geom_bar(stat='identity', position='dodge', width=0.003,alpha=1.0, fill='black') +
  labs(title = "Py", x = "area fraction", y = "TWI")+
  coord_cartesian(ylim=c(0,20),xlim=c(0,0.2))

# pre-computed TWI
p3 <- ggplot(df11, aes(x=twi_pre_comp_area, y=twi_pre_comp_value)) +
  geom_bar(stat='identity', position='dodge', width=0.003,alpha=1.0,fill='black') +
  labs(title = "Pre-comp", x = "area fraction", y = "TWI")+
  coord_cartesian(ylim=c(0,20),xlim=c(0,0.2)) 

ggarrange(p1, p2, p3,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)




