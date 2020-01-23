
#Necessary packages.
library(ggplot2)
library(cowplot)
library(dplyr)
library(rfishbase)

############################################
# The egg mass - temp relationship 
############################################

#Note - this is the code used to combine egg masses from Barneche et al. 2018 & Sibly et al. 2018 with temperature data from Fishbase. Temperature data were accessed from Fishbase in March of 2019.
 
#To convert to mm diameter to wet mass, use same conversion factor as Sibly et al. 2018 (assume eggs are spherical and have the density of water)
#barneche <- read.csv("~Barneche et al. 2018 Glob. Ecol. Biogeog..csv", stringsAsFactor=F) %>% group_by(Species) %>% mutate(barneche_mean_eggsize_g = 0.52*(mean(eggSize_mm)/10)^3) %>% group_by(Species) %>% slice(1L) %>% data.frame() %>% select(Species, barneche_mean_eggsize_g)

#sibly<-read.csv("~Sibly et al. 2018 RSOS.csv", stringsAsFactor=F) %>% mutate(sibly_egg_size_g = Offspring.wet.weight..g, sibly_adult_mass = Adult.wet.weight..g)  %>% tidyr::unite_("Species", c("Genus", "Species"), " ") %>% filter(Class != "Elasmobranchii") %>% select(Class,Species, sibly_egg_size_g, sibly_adult_mass)

#This function accesses temperature data from Fishbase:
#fishbase <- popgrowth() %>% data.frame()%>% group_by(SpecCode) %>% mutate(Temp_mean = mean(Temperature, na.rm=T)) %>% slice(1L) %>% ungroup() %>% data.frame() %>% select(SpecCode, Species,  Temp_mean) %>% filter(!is.na(Temp_mean))

#dat <- full_join(sibly, fishbase, by="Species")
#dat <- full_join(dat, barneche)

#dat <- dat %>% mutate(egg_mass_g = ifelse(is.na(sibly_egg_size_g), barneche_mean_eggsize_g, sibly_egg_size_g)) %>% select(Class, Species,  Temp_mean, egg_mass_g) %>% filter(!is.na(egg_mass_g), !is.na(Temp_mean))



#This is the collated fish egg mass-temperature dataset shown in Fig. 2
#Here and elsewhere in this R file, the .csv files called in the read.csv() function can be found in the /data file on github (https://github.com/daan4786/fish-egg-and-temperature/tree/master/data)
d<-read.csv("~collated fish egg mass-temperature data.csv",stringsAsFactor=F) 
head(d)

modd <- lm(log(Egg.mass..g.)~ Preferred.temp, data = d)
summary(modd)




############################################
#The size and temperature dependencies of egg and larval growth and mortality rates.
############################################


#egg mortality

egg_mort <- read.csv("~Pepin 1991 Can. J. Fish. Aquat. Sci. egg mortality.csv")
head(egg_mort)

egg_mort_mod <- lm(log(egg.mortality.rate..1.day.)~Temp.C, data = egg_mort)
summary(egg_mort_mod)
confint(egg_mort_mod)

be <- exp(egg_mort_mod$coef[1])
Eze <- egg_mort_mod$coef[2]

be_low <- exp(confint(egg_mort_mod)[1,1])
Eze_low <- confint(egg_mort_mod)[2,1]

be_high <- exp(confint(egg_mort_mod)[1,2])
Eze_high <- confint(egg_mort_mod)[2,2]

#egg development

egg_dev<-read.csv("~Pauly & Pullin 1988 Enviro. Biol. Fish. egg development time.csv", header=T)
str(egg_dev)

egg_dev_mod <- lm(ln.dt~ln.egg.mass.g+temp.c, data = egg_dev)
summary(egg_dev_mod)

ae <- exp(egg_dev_mod$coef[1])
Ege <- egg_dev_mod$coef[3]
ne <- egg_dev_mod$coef[2]

ae_low <- exp(confint(egg_dev_mod)[1,1])
Ege_low <- confint(egg_dev_mod)[3,1]
ne_low <- confint(egg_dev_mod)[2,1]

ae_high <- exp(confint(egg_dev_mod)[1,2])
Ege_high <- confint(egg_dev_mod)[3,2]
ne_high <- confint(egg_dev_mod)[2,2]

#larval growth

larv_growth <- read.csv("~houde 1989 Bull. Mar. Sci. larval growth.csv", stringsAsFactor = T)
head(larv_growth)

l_growth_mod <- lm(ln.gr~Temperature+ln.mass, data = larv_growth)
summary(l_growth_mod)

al <- exp(l_growth_mod$coef[1])
Eg <- l_growth_mod$coef[2]
nl <- l_growth_mod$coef[3]

al_low <- exp(confint(l_growth_mod)[1,1])
Eg_low <- confint(l_growth_mod)[2,1]
nl_low <- confint(l_growth_mod)[3,1]

al_high <- exp(confint(l_growth_mod)[1,2])
Eg_high <- confint(l_growth_mod)[2,2]
nl_high <- confint(l_growth_mod)[3,2]


#larval mortality

larv_mort <- read.csv("~Pepin 1991 Can. J. Fish. Aquat. Sci. larval mortality.csv", stringsAsFactor=F)
str(larv_mort)

l_mort_mod <- lm(log(mort.rate.1.day)~ln.mass.g+Temp, data = larv_mort)
summary(l_mort_mod)

bl <- exp(l_mort_mod$coef[1])
Ez <- l_mort_mod$coef[3]
xl <- l_mort_mod$coef[2]

bl_low <- exp(confint(l_mort_mod)[1,1])
Ez_low <- confint(l_mort_mod)[3,1]
xl_low <- confint(l_mort_mod)[2,1]

bl_high <- exp(confint(l_mort_mod)[1,2])
Ez_high <- confint(l_mort_mod)[3,2]
xl_high <- confint(l_mort_mod)[2,2]

#Mass at hatch

larv_hat <- read.csv("~/Documents/research/working projects/dimensionless life history numbers/working code, data and writing/off size/temperature/table 3 Pepin 1991 egg diam-LHAT.csv", stringsAsFactor=F)
head(larv_hat)

l_hat <- lm(log(mass.at.hatch)~log(egg.mass), data = larv_hat)
summary(l_hat)

I<-exp(l_hat$coef[1])
q<-(l_hat$coef[2])

I_low<-exp(confint(l_hat)[1,1])
q_low<-confint(l_hat)[2,1]

I_high<-exp(confint(l_hat)[1,2])
q_high<-confint(l_hat)[2,2]

############################################
#The model
############################################

my_mod <- function(m, Te){
	
	#Using fitted value for each parameter:
	survival_best <- exp( ( (bl/al) * (1/(xl-nl +1)) *exp((Ez-Eg)*Te)* ((I*m^q)^(xl-nl+1) - (0.01)^(xl-nl+1)))  - ((be * ae) * exp((Eze+Ege)*Te)*m^(ne)))  * (1/m)
	
	#Using lower 95%CI for each parameter:
	survival_low <- exp( ( (bl_low/al_low) * (1/(xl_low-nl_low +1)) *exp((Ez_low-Eg_low)*Te) * ((I_low*m^q_low)^(xl_low-nl_low+1) - (0.01)^(xl_low-nl_low+1)))  - ((be_low * ae_low) * exp((Eze_low+Ege_low)*Te)*m^(ne_low)))  * (1/m)
	
	#Using higher 95%CI for each parameter:
	survival_high <- exp( ( (bl_high/al_high) * (1/(xl_high-nl_high +1))*exp((Ez_high-Eg_high)*Te) * ((I_high*m^q_high)^(xl_high-nl_high+1) - (0.01)^(xl_high-nl_high+1)))  - ((be_high * ae_high) * exp((Eze_high+Ege_high)*Te)*m^(ne_high)))  * (1/m)
	
	d<-data.frame(cbind(survival_best, survival_low, survival_high))
	
	return(d)
	
}


############################################
#Make Fig. 2
############################################


mass <- seq(0.000001, 10, by = 0.000001)

temp <- seq(0, 30, by = 1)

temp_dep <- data.frame(temperature = temp, optimal_mass_best = rep(0, length(temp)), optimal_mass_low = rep(0, length(temp)), optimal_mass_high = rep(0, length(temp)), survival = rep(0, length(temp)))

#For each temperature, I find the egg mass that maximizes survivorship per clutch mass, described by the my_mod() function. This egg mass is the optimal mass, and is stored in the temp_dep df. 
for(i in 1:length(temp)){
	
	d_d_best <- data.frame(l= mass, s= my_mod(mass,   temp[i])[,1], ss = (my_mod(mass,   temp[i])[,1]*mass))
	d_d_low <- data.frame(l= mass, s= my_mod(mass,   temp[i])[,2])
	d_d_high <- data.frame(l= mass, s= my_mod(mass,   temp[i])[,3])

	temp_dep$optimal_mass_best[i] <- d_d_best$l[d_d_best$s == max(d_d_best$s )]
	temp_dep$optimal_mass_low[i] <- d_d_low$l[d_d_low$s == max(d_d_low$s )]
	temp_dep$optimal_mass_high[i] <- d_d_high$l[d_d_high$s == max(d_d_high$s )]
	temp_dep$survival[i] <- d_d_best$ss[d_d_best$s == max(d_d_best$s )]

}

#This is the predicted relationship
mod <- lm(log(optimal_mass_best) ~ temperature, data = temp_dep)
summary(mod)

#Plot the observed relationship from the "d" df and the predicted relationship on top of one another

d %>% ggplot(aes(x = Temp_mean, y = log(egg_mass_g) )) + geom_point( shape = 21, stroke=1.1) + geom_smooth(method = "lm", se = F, color = "black") +  ylab("ln egg mass (g)") + xlab(expression(paste("temperature (", ""^{o}, "C)")))  + theme(legend.position = c(0.5, 0.9), legend.title = element_blank()) + ylim(c(-12, -1))   + geom_line(data = temp_dep, aes(x = temperature, y = log(optimal_mass_best)), size = 1,  color = "red")   + geom_line(data = temp_dep, aes(x = temperature, y = log(optimal_mass_low)), size = 0.5, linetype=2,  color = "red")  + geom_line(data = temp_dep, aes(x = temperature, y = log(optimal_mass_high)), size = 0.5, linetype=2,  color = "red") + annotate("text", x = 6, y = -11.5, label = "y == 0.0029 * e^{-0.09*x}", parse=T, size = 5) + annotate("text", x = 6, y = -10.5, label = "y == 0.011 * e^{-0.09*x}", parse=T, size = 5, color = "red")


############################################
#Make Fig. 1
############################################

mass <- seq(0.00001, 0.1, by = 0.00001) 

#Illustrate how the fitness-egg size curves change with temperatures (0, 5, 10, & 15 C)
surv_data <- data.frame(mass = c(mass,mass,mass,mass), temp = c(rep(0, length(mass)), rep(5, length(mass)), rep(10, length(mass)), rep(15, length(mass))) ) %>% rowwise() %>% mutate(survival_per_fecund_best = my_mod(mass,  temp)[,1], survival_per_fecund_low = my_mod(mass,  temp)[,2],survival_per_fecund_high = my_mod(mass,  temp)[,3]) %>% data.frame()
head(surv_data)

#The plot
surv_data %>% ggplot(aes(x = log(mass), y = log(survival_per_fecund_best))) + geom_line(aes(color = as.factor(temp)), size = 1)  + theme( legend.position = c(0, 0.9)) + labs(color = expression(paste("Temp ", ""^{o}, "C"))) + guides(color = guide_legend(ncol=2)) + scale_color_manual(values = c("blue", "green2", "red", "orange")) + xlab(expression(paste("ln"," egg mass (g)")))+ ylab(expression(paste("ln"," number surviving per clutch mass"))) + geom_point(x = log(temp_dep$optimal_mass_best[1]), y = log(max(filter(surv_data, temp==0)$survival_per_fecund_best)), color = "blue", size = 2) + geom_point(x = log(temp_dep$optimal_mass_best[6]), y = log(max(filter(surv_data, temp==5)$survival_per_fecund_best)), color = "green2", size = 2) + geom_point(x = log(temp_dep$optimal_mass_best[11]), y = log(max(filter(surv_data, temp==10)$survival_per_fecund_best)), color = "red", size = 2) + geom_point(x = log(temp_dep$optimal_mass_best[16]), y = log(max(filter(surv_data, temp==15)$survival_per_fecund_best)), color = "orange", size = 2) + xlim(c(-9,-2)) +ylim(c(-2, 2))



############################################
#Make Fig. 3
############################################

#Plot size and temperature dependence of egg and larval mortality and growth rates at two temperatures, to help visualize why smaller eggs are favored in warmer environments. 

rates <- function(m, t){
	
	mass <- m
	temp <- t
	dt <- 1/(65*m^0.2*exp(-0.09*t)) 
	emr <- 0.032*exp(0.17*t)
	lgr <- 0.048*exp(0.063*t)
	lmr <- 0.014*m^-0.23*exp(0.063*t)	
	return(cbind(mass, temp, dt, emr, lgr, lmr))
	
}
mass <- seq(0.0001, 0.1, by = 0.0001)

#Make a df to help plot.
rd <- data.frame(mass = c(mass,mass), temp = c(rep(5, length(mass)), rep(10, length(mass))) ) %>% rowwise() %>% mutate(dt = rates(mass,  temp)[,3], emr = rates(mass,  temp)[,4], lgr = rates(mass,  temp)[,5], lmr = rates(mass,  temp)[,6]) %>% data.frame()

#Make plots.
larvae_hot<-
ggplot() + geom_line(data = filter(rd, temp == 10), aes(x = mass, y = lmr, linetype="1"), color = "red", size = 1)  + geom_line(data = filter(rd, temp == 10), aes(x = mass, y = lgr, linetype="2"), color = "red",  size = 1)   +  ggtitle(expression(paste("larvae at ", 10^{o}, "C")))   + geom_vline(xintercept = temp_dep$optimal_mass_best[temp_dep$temperature == 10],  linetype = "dotted", size = 0.8, color = "red") +  ylab("") + xlab("")  + scale_linetype_manual(name = "",values = c("1"=1,"2"=2), labels = c("mortality", "growth"))  + theme(legend.position=c(0.5,0.95), legend.key.width = unit(1, "cm"), plot.title = element_text(face = "plain")) + scale_x_continuous(breaks = c(0, 0.005, 0.01), limits = c(0, 0.01)) + scale_y_continuous(breaks = c(0.05, 0.15, 0.25), limits = c(0.05, 0.25))

egg_hot<-
ggplot() + geom_line(data = filter(rd, temp == 10), aes(x = mass, y = dt, linetype="2"), color = "red", size = 1) + geom_line(data = filter(rd, temp == 10), aes(x = mass, y = emr, linetype="1"), color = "red", size = 1)   +  ggtitle(expression(paste("eggs at ", 10^{o}, "C")))   + geom_vline(xintercept = temp_dep$optimal_mass_best[temp_dep$temperature == 10], color = "red", linetype = "dotted", size = 0.8)   + ylab("") + xlab("") + scale_linetype_manual(name = "",values = c("1"=1,"2"=2), labels = c("mortality", "development"))  + theme(legend.position=c(0.5,0.95), legend.key.width = unit(1, "cm"), plot.title = element_text(face = "plain")) + scale_x_continuous(breaks = c(0, 0.005, 0.01), limits = c(0, 0.01))+ scale_y_continuous(breaks = c(0.05, 0.15, 0.25), limits = c(0.05, 0.25))

larvae<-
ggplot() + geom_line(data = filter(rd, temp == 5), aes(x = mass, y = lmr), color = "blue", size = 1)  + geom_line(data = filter(rd, temp == 5), aes(x = mass, y = lgr), color = "blue", linetype=2, size = 1)  + theme(legend.position="none", plot.title = element_text(face = "plain"))   +  ggtitle(expression(paste("larvae at ", 5^{o}, "C")))   + geom_vline(xintercept = temp_dep$optimal_mass_best[temp_dep$temperature == 5], color = "blue", linetype = "dotted", size = 0.8)  +  ylab("") + xlab("") + scale_x_continuous(breaks = c(0, 0.005, 0.01), limits = c(0, 0.01))+ scale_y_continuous(breaks = c(0.05, 0.15, 0.25), limits = c(0.05, 0.25))

egg<-
ggplot() + geom_line(data = filter(rd, temp == 5), aes(x = mass, y = dt), color = "blue", linetype=2, size = 1) + geom_line(data = filter(rd, temp == 5), aes(x = mass, y = emr), color = "blue", size = 1) + scale_color_manual(values = c("green2", "orange")) + theme(legend.position="none", plot.title = element_text(face = "plain"))  +  ggtitle(expression(paste("eggs at ", 5^{o}, "C")))  + geom_vline(xintercept = temp_dep$optimal_mass_best[temp_dep$temperature == 5], color = "blue", linetype = "dotted", size = 0.8)   + ylab("") + xlab("") + scale_x_continuous(breaks = c(0, 0.005, 0.01), limits = c(0, 0.01))+ scale_y_continuous(breaks = c(0.05, 0.15, 0.25), limits = c(0.05, 0.25))

#Put plots together.
p<-plot_grid(egg_hot, larvae_hot, egg, larvae, align =c("l", "r" ), labels = "AUTO", label_x=0.23, label_y=0.89, label_fontface="plain")
ggdraw(p ) + draw_label("mass (g)", x = 0.5, y = 0.03, size = 14) + draw_label("rate (1/day)", x = 0.03, y = 0.5, angle = 90, size = 14)




