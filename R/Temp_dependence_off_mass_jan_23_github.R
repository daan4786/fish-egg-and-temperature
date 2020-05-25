
#Necessary packages.
library(ggplot2)
library(cowplot)
library(dplyr)
library(rfishbase)
library(MCMCglmm)
library(HDInterval)

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


egg_mort_mod_mcmc <- MCMCglmm(log(egg.mortality.rate..1.day.)~Temp.C, data = egg_mort, nitt=100000, burnin=10000)
summary(egg_mort_mod_mcmc)
#plot(egg_mort_mod_mcmc) #inspect trace and density plots to ensure convergence.


posterior_egg_mort_temp_dep <- egg_mort_mod_mcmc$Sol[,2]
posterior_egg_mort_int <- egg_mort_mod_mcmc$Sol[,1]


be <- exp(mean(posterior_egg_mort_int))
Eze <- mean(posterior_egg_mort_temp_dep)


#egg development

egg_dev<-read.csv("~Pauly & Pullin 1988 Enviro. Biol. Fish. egg development time.csv", header=T)


egg_dev_mod_mcmc <- MCMCglmm(ln.dt~ln.egg.mass.g+temp.c, data = egg_dev, nitt=100000, burnin=10000)
summary(egg_dev_mod_mcmc)
#plot(egg_dev_mod_mcmc) #inspect trace and density plots to ensure convergence.


posterior_egg_dev_mass_dep <- egg_dev_mod_mcmc$Sol[,2]
posterior_egg_dev_temp_dep <- egg_dev_mod_mcmc$Sol[,3]
posterior_egg_dev_int <- egg_dev_mod_mcmc$Sol[,1]

ae <- exp(mean(posterior_egg_dev_int))
Ege <- mean(posterior_egg_dev_temp_dep)
ne <- mean(posterior_egg_dev_mass_dep)

#larval growth

larv_growth <- read.csv("~houde 1989 Bull. Mar. Sci. larval growth.csv", stringsAsFactor = T)
head(larv_growth)

l_growth_mod_mcmc <- MCMCglmm(ln.gr~temp+ln.mass, data = barn_l_g, nitt=100000, burnin=10000)
summary(l_growth_mod_mcmc)
#plot(l_growth_mod_mcmc)

posterior_l_growth_mass_dep <- l_growth_mod_mcmc$Sol[,3]
posterior_l_growth_temp_dep <- l_growth_mod_mcmc$Sol[,2]
posterior_l_growth_int <- l_growth_mod_mcmc$Sol[,1]

al <- exp(mean(posterior_l_growth_int))
Eg <- mean(posterior_l_growth_temp_dep)
nl <- mean(posterior_l_growth_mass_dep)


#larval mortality

larv_mort <- read.csv("~Pepin 1991 Can. J. Fish. Aquat. Sci. larval mortality.csv", stringsAsFactor=F)


l_mort_mod_mcmc <- MCMCglmm(log(mort.rate.1.day)~ln.mass.g+Temp, data = larv_mort, nitt=100000, burnin=10000)
summary(l_mort_mod_mcmc)
#plot(l_mort_mod_mcmc)

posterior_l_mort_mass_dep <- l_mort_mod_mcmc$Sol[,2]
posterior_l_mort_temp_dep <- l_mort_mod_mcmc$Sol[,3]
posterior_l_mort_int <- l_mort_mod_mcmc$Sol[,1]


bl <- exp(mean(posterior_l_mort_int))
Ez <- mean(posterior_l_mort_temp_dep)
xl <- mean(posterior_l_mort_mass_dep)

#Mass at hatch

larv_hat <- read.csv("~Pepin 1991 Can. J. Fish. Aquat. Schi. egg mass-mass at hatch.csv", stringsAsFactor=F)

l_met_mcmc <- MCMCglmm(log(mass.at.hatch)~log(egg.mass), data = larv_met, nitt=100000, burnin=10000)
summary(l_met_mcmc)
#plot(l_met_mcmc)

posterior_l_met_temp_dep <- l_met_mcmc$Sol[,2]
posterior_l_met_int <- l_met_mcmc$Sol[,1]

I<-exp(mean(posterior_l_met_int))
q<-mean(posterior_l_met_temp_dep)


############################################
#The model
############################################

my_mod <- function(m, Te, bl, al, xl, nl, I, q, be, ae, Eze, Ege, ne, Ez, Eg){
	
	survival_best <- exp( ( (bl/al) * (1/(xl-nl +1)) * exp((Ez - Eg)*Te) * ((I*m^q)^(xl-nl+1) - (0.1)^(xl-nl+1)))  - ((be * ae) * exp((Eze+Ege)*Te)*m^(ne)))  * (1/m)
	
	d<-data.frame(survival_best)
	
	return(d)
	
}


############################################
#Make Fig. 2
############################################


mass000000 <- seq(0.00000000000000000000000000001, 0.00000000000000000000000001, by = 0.00000000000000000000000000001)
mass00000 <- seq(0.00000000000000000000000001, 0.00000000000000000000001, by = 0.00000000000000000000000001)
mass0000 <- seq(0.00000000000000000000001, 0.00000000000000000001, by = 0.00000000000000000000001)
mass000 <- seq(0.00000000000000000001, 0.00000000000000001, by = 0.00000000000000000001)
mass00 <- seq(0.00000000000000001, 0.00000000000001, by = 0.00000000000000001)
mass0 <- seq(0.00000000000001, 0.00000000001, by = 0.00000000000001)
mass1 <- seq(0.00000000001, 0.00000001, by = 0.00000000001)
mass2 <- seq(0.00000001, 0.00001, by = 0.00000001)
mass3 <- seq(0.00001, 0.01, by = 0.00001)
mass4 <- seq(0.01, 10, by = 0.01)
mass5 <- seq(10, 100, by = 1)
mass <- c(mass000000, mass00000, mass0000, mass00, mass0, mass1, mass2, mass3, mass4, mass5)

temp <- seq(0, 30, by = 1)

temp_dep <- data.frame(temperature = temp, optimal_mass_best = rep(0, length(temp)), optimal_mass_low = rep(0, length(temp)), optimal_mass_high = rep(0, length(temp)), survival = rep(0, length(temp)))


size_pd <- 50000

posterior_distribution_predictions <- data.frame(temp = c(rep(0, size_pd), rep(1, size_pd),rep(2, size_pd),rep(3,(size_pd)),rep(4,(size_pd)),rep(5,(size_pd)),rep(6,(size_pd)),rep(7,(size_pd)),rep(8,(size_pd)),rep(9,(size_pd)),rep(10,(size_pd)),rep(11,(size_pd)),rep(12,(size_pd)),rep(13,(size_pd)),rep(14,(size_pd)),rep(15,(size_pd)),rep(16,(size_pd)),rep(17,(size_pd)),rep(18,(size_pd)),rep(19,(size_pd)),rep(20,(size_pd)),rep(21,(size_pd)),rep(22,(size_pd)),rep(23,(size_pd)),rep(24,(size_pd)),rep(25,(size_pd)),rep(26,(size_pd)),rep(27,(size_pd)),rep(28,(size_pd)),rep(29,(size_pd)),rep(30,(size_pd))), optimal_mass_pred = rep(0, (size_pd*length(temp))))

#For each temperature, I find the egg mass that maximizes survivorship per clutch mass, described by the my_mod() function. This egg mass is the optimal mass, and is stored in the temp_dep df. 
for(i in 1:length(temp)){
	
	#Generate a df of parental fitness ("s") vs. mass ("l"). 
	d_d_best <- data.frame(l= mass, s= my_mod(mass, temp[i],bl, al, xl, nl, I, q, be, ae, Eze, Ege, ne, Ez, Eg)[,1], ss = (my_mod(mass,temp[i],bl, al, xl, nl, I, q, be, ae, Eze, Ege, ne, Ez, Eg)[,1]*mass))
	
	#Find and store the mass ("l") that maximizes parental fitness ("s").
	temp_dep$optimal_mass_best[i] <- d_d_best$l[d_d_best$s == max(d_d_best$s )]
	temp_dep$survival[i] <- d_d_best$ss[d_d_best$s == max(d_d_best$s )]
	
	x <- rep(0, size_pd)
	
	for(j in 1:size_pd){
	
		#At each temperature, I generate "size_pd" (=50,000 in this case) predictions for optimal egg size by randomly sampling the posterior distribution for each model parameter.	
		d_d_pd <- data.frame(l= mass, s= my_mod(mass, temp[i], exp(sample(posterior_l_mort_int, 1)), exp(sample(posterior_l_growth_int, 1)), sample(posterior_l_mort_mass_dep, 1), sample(posterior_l_growth_mass_dep, 1), exp(sample(posterior_l_met_int, 1)), sample(posterior_l_met_temp_dep, 1), exp(sample(posterior_egg_mort_int, 1)), exp(sample(posterior_egg_dev_int, 1)), sample(posterior_egg_mort_temp_dep, 1), sample(posterior_egg_dev_temp_dep, 1), sample(posterior_egg_dev_mass_dep, 1), sample(posterior_l_mort_temp_dep, 1), sample(posterior_l_growth_temp_dep, 1))[,1])
	
		x[j] <- d_d_pd$l[d_d_pd$s == max(d_d_pd$s )]
		
	}

posterior_distribution_predictions$optimal_mass_pred[posterior_distribution_predictions$temp == (i-1)] <- x	

print(i)
	
}

#This is the predicted relationship
mod <- lm(log(optimal_mass_best) ~ temperature, data = temp_dep)
summary(mod)


#This df gives the lower and upper 95% Highest Density Intervals for the 50,000 predictions which were generated by randomly sampling the posterior distributions for each model parameter.
uncertainty_df <- posterior_distribution_predictions %>% group_by(temp) %>% summarize( ci_low = hdi(log(optimal_mass_pred))[1], ci_high = hdi(log(optimal_mass_pred))[2])

#Plot the observed relationship from the "d" df and the predicted relationship on top of one another. Also include the lower and upper 95% HDIs from the uncertainty df. 
d %>% ggplot(aes(x = Preferred.temp, y = log(Egg.mass..g.) )) + geom_point( shape = 21, stroke=1.1) + geom_smooth(method = "lm", se = F, color = "black") +  ylab("ln egg mass (g)") + xlab(expression(paste("temperature (", ""^{o},"C)")))  + theme(legend.position = c(0.5, 0.9), legend.title = element_blank())    + geom_line(data = temp_dep, aes(x = temperature, y = log(optimal_mass_best)), size = 1,  color = "red")   + geom_line(data = uncertainty_df, aes(x = temp, y = ci_low), method = "lm", color = "red", linetype = 2)  + geom_line(data = uncertainty_df, aes(x = temp, y = ci_high), color = "red", linetype = 2) + annotate("text", x = 4.5, y = -23.5, label = "y == 0.003 * e^{-0.09*x}", parse=T, size = 5) + annotate("text", x = 4.5, y = -22, label = "y == 0.013 * e^{-0.11*x}", parse=T, size = 5, color = "red") 


############################################
#Make Fig. 1
############################################


mass2 <- seq(0.00000001, 0.00001, by = 0.00000001)
mass3 <- seq(0.00001, 0.01, by = 0.00001)
mass4 <- seq(0.01, 10, by = 0.01)
mass5 <- seq(10, 100, by = 1)
mass <- c(mass3, mass4, mass5)


#Illustrate how the fitness-egg size curves change with temperatures (0, 5, 10, & 15 C)
surv_data <- data.frame(mass = c(mass,mass,mass,mass), temp = c(rep(0, length(mass)), rep(5, length(mass)), rep(10, length(mass)), rep(15, length(mass))) ) %>% rowwise() %>% mutate(survival_per_fecund_best = my_mod(mass,  temp, bl, al, xl, nl, I, q, be, ae, Eze, Ege, ne, Ez, Eg)[,1]) %>% data.frame()
head(surv_data)

#The plot
surv_data %>% ggplot(aes(x = log(mass), y = log(survival_per_fecund_best))) + geom_line(aes(color = as.factor(temp)), size = 1)  + theme( legend.position = c(0, 0.9)) + labs(color = expression(paste("Temp ", ""^{o}, "C"))) + guides(color = guide_legend(ncol=2)) + scale_color_manual(values = c("blue", "green2", "red", "orange")) + xlab("ln egg mass (g)")+ ylab("ln number surviving per clutch mass") + geom_point(x = log(temp_dep$optimal_mass_best[1]), y = log(max(filter(surv_data, temp==0)$survival_per_fecund_best)), color = "blue", size = 2) + geom_point(x = log(temp_dep$optimal_mass_best[6]), y = log(max(filter(surv_data, temp==5)$survival_per_fecund_best)), color = "green2", size = 2) + geom_point(x = log(temp_dep$optimal_mass_best[11]), y = log(max(filter(surv_data, temp==10)$survival_per_fecund_best)), color = "red", size = 2) + geom_point(x = log(temp_dep$optimal_mass_best[16]), y = log(max(filter(surv_data, temp==15)$survival_per_fecund_best)), color = "orange", size = 2) + xlim(c(-10, 0)) + ylim(c(-6, 0))



############################################
#Make Fig. 3
############################################

#Plot size and temperature dependence of egg and larval mortality and growth rates at two temperatures, to help visualize why smaller eggs are favored in warmer environments. 


rates <- function(m, t){
	
	mass <- m
	temp <- t
	dt <- 1/(ae*m^ne*exp(Ege*t))
	emr <- be*exp(Eze*t)
	lgr <- al*(m^nl)*exp(Eg*t)/m
	lmr <- bl*(m^xl)*exp(Ez*t)	
	return(cbind(mass, temp, dt, emr, lgr, lmr))
	
}
mass <- seq(0.0001, 0.1, by = 0.0001)


#Make a df to help plot.
rd <- data.frame(mass = c(mass,mass), temp = c(rep(5, length(mass)), rep(10, length(mass))) ) %>% rowwise() %>% mutate(dt = rates(mass,  temp)[,3], emr = rates(mass,  temp)[,4], lgr = rates(mass,  temp)[,5], lmr = rates(mass,  temp)[,6]) %>% data.frame()

#Make plots.

larvae_hot<-
ggplot() + geom_line(data = filter(rd, temp == 10), aes(x = mass, y = lmr, linetype="1"), color = "red", size = 1)  + geom_line(data = filter(rd, temp == 10), aes(x = mass, y = lgr, linetype="2"), color = "red",  size = 1)   +  ggtitle(expression(paste("larvae at ", 10^{o}, "C")))   + geom_vline(xintercept = temp_dep$optimal_mass_best[temp_dep$temperature == 10],  linetype = "dotted", size = 0.8, color = "red") +  ylab("") + xlab("")  + scale_linetype_manual(name = "",values = c("1"=1,"2"=2), labels = c("mortality", "growth"))  + theme(legend.position=c(0.5,0.95), legend.key.width = unit(1, "cm"), plot.title = element_text(face = "plain")) + scale_x_continuous(breaks = c(0,  0.01, 0.02), limits = c(0, 0.02)) + scale_y_continuous(breaks = c(0.05, 0.15, 0.25), limits = c(0.05, 0.25))

egg_hot<-
ggplot() + geom_line(data = filter(rd, temp == 10), aes(x = mass, y = dt, linetype="2"), color = "red", size = 1) + geom_line(data = filter(rd, temp == 10), aes(x = mass, y = emr, linetype="1"), color = "red", size = 1)   +  ggtitle(expression(paste("eggs at ", 10^{o}, "C")))   + geom_vline(xintercept = temp_dep$optimal_mass_best[temp_dep$temperature == 10], color = "red", linetype = "dotted", size = 0.8)   + ylab("") + xlab("") + scale_linetype_manual(name = "",values = c("1"=1,"2"=2), labels = c("mortality", "development"))  + theme(legend.position=c(0.5,0.95), legend.key.width = unit(1, "cm"), plot.title = element_text(face = "plain")) + scale_x_continuous(breaks = c(0,  0.01, 0.02), limits = c(0, 0.02))+ scale_y_continuous(breaks = c(0.05, 0.15, 0.25), limits = c(0.05, 0.25))

larvae<-
ggplot() + geom_line(data = filter(rd, temp == 5), aes(x = mass, y = lmr), color = "blue", size = 1)  + geom_line(data = filter(rd, temp == 5), aes(x = mass, y = lgr), color = "blue", linetype=2, size = 1)  + theme(legend.position="none", plot.title = element_text(face = "plain"))   +  ggtitle(expression(paste("larvae at ", 5^{o}, "C")))   + geom_vline(xintercept = temp_dep$optimal_mass_best[temp_dep$temperature == 5], color = "blue", linetype = "dotted", size = 0.8)  +  ylab("") + xlab("") + scale_x_continuous(breaks = c(0,  0.01, 0.02), limits = c(0, 0.02))+ scale_y_continuous(breaks = c(0.05, 0.15, 0.25), limits = c(0.05, 0.25))

egg<-
ggplot() + geom_line(data = filter(rd, temp == 5), aes(x = mass, y = dt), color = "blue", linetype=2, size = 1) + geom_line(data = filter(rd, temp == 5), aes(x = mass, y = emr), color = "blue", size = 1) + scale_color_manual(values = c("green2", "orange")) + theme(legend.position="none", plot.title = element_text(face = "plain"))  +  ggtitle(expression(paste("eggs at ", 5^{o}, "C")))  + geom_vline(xintercept = temp_dep$optimal_mass_best[temp_dep$temperature == 5], color = "blue", linetype = "dotted", size = 0.8)   + ylab("") + xlab("") + scale_x_continuous(breaks = c(0,  0.01, 0.02), limits = c(0, 0.02))+ scale_y_continuous(breaks = c(0.05, 0.15, 0.25), limits = c(0.05, 0.25))

#put plots together.
p<-plot_grid(egg_hot, larvae_hot, egg, larvae, align =c("l", "r" ), labels = "AUTO", label_x=0.23, label_y=0.89, label_fontface="plain")
ggdraw(p ) + draw_label("mass (g)", x = 0.5, y = 0.03, size = 14) + draw_label("rate (1/day)", x = 0.03, y = 0.5, angle = 90, size = 14)


