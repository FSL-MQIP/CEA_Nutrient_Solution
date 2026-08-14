#Challenge study Analysis for CEA project
#Load  packages
library("tidyverse")
library(lmerTest)
library(emmeans)
library(nlsMicrobio)
library(grid)
library(gridExtra)
library(extrafont)
font_import() # Make sure to type 'Y' when prompted
loadfonts(device = "win", quiet = TRUE)

#save the R environment
save.image(file='C:/Users/tmwal/OneDrive/Desktop/Challengestudydata04242026.RData')

#See the following link for help with interpreting linear models, regression and anova
#https://cscu.cornell.edu/workshop/interpreting-linear-models-regression-and-anova/

#load saved R environment
load("C:/Users/tmwal/OneDrive/Desktop/Challengestudydata04242026.RData")

#import Total Plate Count (TPC) files 
TPC<-read.csv('/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Induvidual_Water_Samples_for_manuscipt_07022025.csv')
TPC<-TPC%>%na.omit()
TPC$Pond_Type_Number <- paste(TPC$Pond_Type,TPC$Pond_Number)

#calculate TPC average values
TPC1 <- TPC %>% group_by(Sample_Name) %>% mutate(Average_Conc=mean(Concentration))
TPC1 <- TPC1 %>% group_by(Sample_Name) %>% mutate(Log_Average_Conc=log10(Average_Conc))%>%ungroup

#Use lmer to make mixed model to test impact of pond type on log count data using unavgd replicates
  # Pond Number:fixed effect
  # Location:fixed effect
TPC2<-TPC1%>%distinct(Pond_Type_Number,Pond_Type,Pond_Number, Log_Average_Conc, .keep_all = TRUE)
TPC2$Pond_Type_Number <- as.factor(TPC2$Pond_Type_Number)

TPC_lmer <- lme4::lmer(log10(Average_Conc) ~ Pond_Type+ Location+(1|Pond_Type_Number), TPC2)
summary(TPC_lmer)
#for pond type
emmeans(TPC_lmer, pairwise ~ Pond_Type, type = "response")
#To calculate % change on your own, Subtract 1 from the factor value (the ratio column in the 
#emmeans output), then multiply by 100

#Conventional / Organic
(0.354-1)*100 #-64.6% change; decr by ~65%


#for location
emmeans(TPC_lmer, pairwise ~  Location, type = "response")

#location H / M
(1.28-1)*100 #28% change; incr by 28%
#location H / S
(1.39-1)*100 #39% change; incr by ~39%
#location M / S
(1.08-1)*100 #8% change; incr by ~8%

# Simple Linear Model
TPC_final <- lm(Log_Average_Conc ~ Pond_Type + Location, data = TPC2)

# Check diagnostic plots 
plot(TPC_final)
anova(TPC_final)
summary(TPC_final)

#Calculate the average log concentration of conventional vs organic samples
TPC3 <- TPC2 %>% group_by(Pond_Type) %>% 
  mutate(Average_Conc_byPondType=mean(Average_Conc))%>%
  mutate(Log_Average_Conc_byPondType=log10(Average_Conc_byPondType))%>%
  mutate(log_SD_Average_Conc_byPond=sd(Log_Average_Conc)) %>%ungroup

TPC_logavgconcbypondtype<-TPC3%>%distinct(Pond_Type, Log_Average_Conc_byPondType, .keep_all = TRUE)


#plot

apccolor <- rep(c( "#D55E00", 
                     "#0072B2"))


#plot
apcplot <- ggplot(TPC3, aes(x = Pond_Type, y = Log_Average_Conc, 
                              color = Pond_Type)) + 
  geom_boxplot(size = 5,width = 0.5) + 
  geom_point(position=position_dodge(width=0.75),color="black", size=10)+
  labs(title = expression( ~ italic("Aerobic Plate Count")),x = "Pond Type", 
       y = expression(Concentration~ (log [10] ~CFU/mL)))+
  theme_bw() +
  scale_color_manual(values = apccolor) +
  scale_shape_manual(values = c(16, 17, 15, 18)) + 
  guides(color = guide_legend(title = "Pond Type"))+
  theme(text = element_text(family = "Times New Roman"),
        strip.text = element_text(size = 40, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
       legend.title = element_text(size = 50), 
        legend.text = element_text(size = 50),
        axis.text.x = element_text(angle = 0,  size = 50),
        axis.text.y = element_text(angle = 0, vjust = 0.3, hjust = 0.5, size = 50, margin = margin(0,0,0,20)),
        axis.title.y = element_text(size = 50),
        axis.title.x = element_text(size = 50),
        plot.title = element_text(size = 55, face = "italic", hjust = 0.5))

apcplot

ggsave(
  filename = "/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/apcplot.tiff", 
  plot = last_plot(),       
  device = "tiff",          
  units = "in", 
  width = 25, 
  height = 15, 
  dpi = 700, 
  compression = "lzw",
  type = "cairo")


#INOCULATION DATA
#import pathogen data
Day0<-read.csv("/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/DAY0_COMBINED_FOR_R_updated_v2.csv")
Day1<-read.csv("/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/DAY1_COMBINED_FOR_R_updated_v2.csv")
Day2<-read.csv("/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/DAY2_COMBINED_FOR_R_updated_v2.csv")
Day5<-read.csv("/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/DAY5_combined_for_R_updated_v2.csv")

#MERGE data all in one frame
df <- rbind(Day0,Day1) 
df <- rbind(df,Day2) 
df <- rbind(df,Day5) 

#add new groupings/ sample name combos
df$Sample_ID_Aliquot<-paste(df$Day,df$Pond_Type_Number,df$Aliquot,sep = '_')
df$Sample_ID<-paste(df$Day,df$Pond_Type_Number,sep = '_')
df$Pond_Type_Number<-as.character(df$Pond_Type_Number)

#Listeria monocytogenes______________________________________________________

#sort out Lmono
Lmono <- df %>% filter(Bacteria=="LMONO")

#Sort out Lmono NC, check that it looks ok
LmonoNC <- df %>% filter(Bacteria=="LMONO")%>% filter(Aliquot=="NC")

#Take out Lmono NC, calc average values
Lmono_A_C <- df %>% filter(Bacteria=="LMONO")%>% filter(Aliquot!="NC")

#Replace -inf in Log_Avg_Conc column with 1.3, which is the Limit of Detection, 
#and replace 0 in Avg_Conc column with 10^1.3
#Lmono_A_C_1$Log_Average_Conc<- gsub("-Inf", "1.3",Lmono_A_C_1$Log_Average_Conc)
x=10^(1.3)
Lmono_A_C$Final_Concentration[Lmono_A_C$Final_Concentration == "0"] <- x
Lmono_A_C$Final_Concentration<- as.numeric(Lmono_A_C$Final_Concentration)

#calculate average values by replicate
Lmono_A_C_byREP <- Lmono_A_C %>% 
  group_by(Day,Pond_Type_Number, Aliquot, Dilution) %>% 
  mutate(Average_Count=mean(Counted_Colonies))%>% 
  mutate(Average_Conc=mean(Final_Concentration)) %>% 
  mutate(Log_Average_Conc=log10(Average_Conc))%>%
  ungroup

Lmono_A_C_1 <- Lmono_A_C_byREP %>%distinct(Day,Pond_Type_Number, Aliquot, .keep_all = TRUE)


#plot

mycolslmono <- rep(c("#E69F00", 
                    "#D55E00", 
                    "#56B4E9",
                    "#0072B2"), 
                  length.out = length(unique(Lmono_A_C_1$Sample_ID)))

Lmono_A_C_1$Day <- factor(Lmono_A_C_1$Day, levels = c("DAY0", "DAY1", "DAY2", "DAY5"), 
                                 labels = c("Day 0", "Day 1", "Day 2","Day 5"))

y.expression <- expression(Concentration~ (log [10] ~CFU/mL))

#Make a continuous day column

Lmono_A_C_1<-Lmono_A_C_1%>%mutate(Day_numeric=Day)%>%mutate(Day_numeric = str_remove_all(Day, "[^0-9.-]")) %>%
  mutate(Day_numeric = as.numeric(Day_numeric))

#plot
r5 <- ggplot(Lmono_A_C_1, aes(x = Day_numeric, y = Log_Average_Conc, 
                                         shape = factor(Pond_Type_Number), 
                                         color = factor(Pond_Type_Number))) + 
  # Map both color and shape to the same factor to merge the legend
  geom_point(size = 18, alpha = 0.9) + 
  labs(title = expression( ~ italic("Listeria monocytogenes")),
       x = "Day", 
       y = y.expression) +  
  theme_bw() +
  # Use show.legend = FALSE to keep the dashed line out of your point legend
  geom_hline(yintercept = 1.3, linetype = "dashed", show.legend = FALSE) +
  ylim(1, 5) +
  # Ensure the manual color scale name matches the shape scale (both are NULL/blank)
  scale_color_manual(values = mycolslmono) +
  # Optional: force shape values if you don't like the defaults
  scale_shape_manual(values = c(16, 17, 15, 18)) + 
  theme(text = element_text(family = "Times New Roman"),
        strip.text = element_text(size = 40, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank(), # This helps merging
        legend.text = element_text(size = 30),
        axis.text.x = element_text(angle = 0, vjust = 0.3, hjust = 1, size = 40, margin = margin(0,0,20,0)),
        axis.text.y = element_text(angle = 0, vjust = 0.3, hjust = 0.5, size = 40, margin = margin(0,0,0,20)),
        axis.title.y = element_text(size = 50),
        axis.title.x = element_text(size = 50),
        plot.title = element_text(size = 55, face = "italic", hjust = 0.5))

r5

ggsave(
  filename = "/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Lmonolin_contxaxis_poster3.tiff", 
  plot = last_plot(),       
  device = "tiff",          
  units = "in", 
  width = 25, 
  height = 15, 
  dpi = 700, 
  compression = "lzw",
  type = "cairo"            
)

#COMPARING LINEAR AND MAFART WEIBULL MODELS
# Fit a linear model
Lmonomodel <- lm(Log_Average_Conc ~ Day_numeric, data = Lmono_A_C_1)

# Generate 4 diagnostic plots (Residuals vs Fitted is first)
plot(Lmonomodel) 
BIC(Lmonomodel)
summary(Lmonomodel)

# Get 95% Confidence Intervals
confint(Lmonomodel)
Lmonomodelparameters_Lm <- summary(Lmonomodel)$coefficients
print(Lmonomodelparameters_Lm)


# To see only the Residuals vs Fitted plot:
plot(Lmonomodel, which = 1)
ggsave("/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Check_linearity/Lmonoresiduals.tiff", units="in", width=20, height=15, dpi=500, compression = 'lzw')

# qqplotplot:
plot(Lmonomodel, which = 2)
ggsave("/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Check_linearity/Lmonoqq.tiff", units="in", width=20, height=15, dpi=500, compression = 'lzw')

# Fit the Lmono data to a weibull distribution (mafart)
mafart_weibull_Lm <- nls(Log_Average_Conc ~ log10_N0 - (Day_numeric/delta)^p, 
                        data =Lmono_A_C_1, 
                        start = list(log10_N0 = 4.661987, delta = 1, p = 0.5))


overview(mafart_weibull_Lm)
plotfit(mafart_weibull_Lm, smooth = TRUE)
BIC(mafart_weibull_Lm)
mafartparameters_Lm <- summary(mafart_weibull_Lm)$coefficients
print(mafartparameters_Lm)

#check residuals
plot(fitted(mafart_weibull_Lm), residuals(mafart_weibull_Lm))
abline(h = 0, lty = 2, col = "red")

#Use lmer to make mixed model to test impact of pond type on log count data using unaveraged aliquot values
# use log(Average_Conc) in the formula to inform future emmeans code that this data is log transformed
# see this vignette for a full explanation: https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html

class(Lmono_A_C_1)
Lmono_lmer <- lmerTest::lmer(log10(Average_Conc) ~ Day * Pond_Type+(1|Pond_Type_Number), Lmono_A_C_1)
summary(Lmono_lmer)
anova(Lmono_lmer)

#add ',type="response" to this emmeans code so it recognizes
#the values are log10 transformed and the responses emmeans provides are exponentiated
#so you dont have to exponentiate by hand
#But you have to calculate % change on your own
emmeans(Lmono_lmer, consec ~ Day,type="response")

#To calculate % change on your own, Subtract 1 from the factor value (the ratio column in the 
#emmeans output), then multiply by 100
#Day1/Day0
(0.6839-1)*100 #-31.61% change; decr by ~32%
#Day2/Day1
(0.0405-1)*100 #-95.95% change; decr by ~96%
#Day2/Day1
(0.0251-1)*100 #-97.49% change; decr by ~97%

#Average the aliquots together to find overall average values by organic/ conventional
Lmono_A_C_bypond <- Lmono_A_C_1%>%
  group_by(Day,Pond_Type) %>% 
  mutate(Average_Log_Conc_byPond=log10(mean(Average_Conc)))%>%
  mutate(log_SD_Average_Conc_byPond=sd(Log_Average_Conc)) %>%
  mutate(log_SER_Average_Conc_byPond=(sd(Log_Average_Conc))/(sqrt(3))) %>%
  ungroup

Lmono_A_C_bypond <- Lmono_A_C_bypond %>%distinct(Day,Pond_Type, .keep_all = TRUE)

Lmono_A_C_bypond_fortable <-Lmono_A_C_bypond%>%select(Sample_ID, Bacteria,Average_Log_Conc_byPond,log_SD_Average_Conc_byPond)%>%
  mutate(Average_Log_Conc_byPond = round(Average_Log_Conc_byPond, digits = 2))%>%
  mutate(log_SD_Average_Conc_byPond = round(log_SD_Average_Conc_byPond, digits = 2))
  
#SALMONELLA________________________________________________________________

#sort out Salm
Salm <- df %>% filter(Bacteria=="SALM")
#sort out Salm NC, check that it looks ok
SalmNC <- df %>% filter(Bacteria=="SALM")%>% filter(Aliquot=="NC")
#take out Salm NC, calc avg values
Salm_A_C <- df %>% filter(Bacteria=="SALM")%>% filter(Aliquot!="NC")

#calculate average values by replicate
Salm_A_C_byREP <- Salm_A_C %>% 
  group_by(Day,Pond_Type_Number, Aliquot, Dilution) %>% 
  mutate(Average_Count=mean(Counted_Colonies))%>% 
  mutate(Average_Conc=mean(Final_Concentration)) %>% 
  mutate(Log_Average_Conc=log10(Average_Conc))%>%
  ungroup

Salm_A_C_1 <- Salm_A_C_byREP %>%distinct(Day,Pond_Type_Number, Aliquot, .keep_all = TRUE)
str(Salm_A_C_1)
Salm_A_C_1$Day<-as.factor(Salm_A_C_1$Day)


#Make a continuous day column

Salm_A_C_1<-Salm_A_C_1%>%mutate(Day_numeric=Day)%>%mutate(Day_numeric = str_remove_all(Day, "[^0-9.-]")) %>%
  mutate(Day_numeric = as.numeric(Day_numeric))

mycolsSalm <- rep(c("#E69F00", 
                    "#D55E00", 
                    "#56B4E9",
                    "#0072B2"), 
                  length.out = length(unique(Salm_A_C_1$Sample_ID)))

Salm_A_C_1$Day <- factor(Salm_A_C_1$Day, levels = c("DAY0", "DAY1", "DAY2", "DAY5"), 
                         labels = c("Day 0", "Day 1", "Day 2","Day 5"))

rsalm <- ggplot(Salm_A_C_1, aes(x = Day_numeric, y = Log_Average_Conc, 
                                         shape = factor(Pond_Type_Number), 
                                         color = factor(Pond_Type_Number))) + 
  # Map both color and shape to the same factor to merge the legend
  geom_point(size = 18, alpha = 0.9) + 
  labs(title = expression(~italic("Salmonella enterica")),
       x = "Day", 
       y = y.expression) +  
  theme_bw() +
  # Use show.legend = FALSE to keep the dashed line out of your point legend
  geom_hline(yintercept = 1.3, linetype = "dashed", show.legend = FALSE) +
  ylim(1, 5) +
  # Ensure the manual color scale name matches the shape scale (both are NULL/blank)
  scale_color_manual(values = mycolsSalm) +
  # Optional: force shape values if you don't like the defaults
  scale_shape_manual(values = c(16, 17, 15, 18)) + 
  theme(text = element_text(family = "Times New Roman"),
        strip.text = element_text(size = 40, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 30),
        axis.text.x = element_text(angle = 0, vjust = 0.3, hjust = 1, size = 40, margin = margin(0,0,20,0)),
        axis.text.y = element_text(angle = 0, vjust = 0.3, hjust = 0.5, size = 40, margin = margin(0,0,0,20)),
        axis.title.y = element_text(size = 50),
        axis.title.x = element_text(size = 50),
        plot.title = element_text(size = 55, face = "italic", hjust = 0.5))

rsalm

ggsave(
  filename = "/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Salml_contxaxis_poster.tiff", 
  plot = last_plot(),     
  device = "tiff",          
  units = "in", 
  width = 25, 
  height = 15, 
  dpi = 700, 
  compression = "lzw",
  type = "cairo")

#COMPARING LINEAR AND MAFART WEIBULL MODELS

# Fit a model
Salmmodel <- lm(Log_Average_Conc ~ Day_numeric, data = Salm_A_C_1)
BIC(Salmmodel)

# To see only the Residuals vs Fitted plot:
plot(Salmmodel, which = 1)
ggsave("/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Check_linearity/Salmresiduals.tiff", units="in", width=20, height=15, dpi=500, compression = 'lzw')

# qqplotplot:
plot(Salmmodel, which = 2)
ggsave("/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Check_linearity/Salmqq.tiff", units="in", width=20, height=15, dpi=500, compression = 'lzw')


#divide Salm_A_C_1 in to organic and conventional
Salm_A_C_1_O<-Salm_A_C_1%>%filter(Pond_Type=="Organic")
Salm_A_C_1_C<-Salm_A_C_1%>%filter(Pond_Type=="Conventional")

#linear model
Salmmodel_C <- lm(Log_Average_Conc ~ Day_numeric, data = Salm_A_C_1_C)
BIC(Salmmodel_C)
Salmmodel_O <- lm(Log_Average_Conc ~ Day_numeric, data = Salm_A_C_1_O)
BIC(Salmmodel_O)

#mafart conventional
mafart_weibull_SC <- nls(Log_Average_Conc ~ log10_N0 - (Day_numeric/delta)^p, 
                      data =Salm_A_C_1_C, 
                      start = list(log10_N0 = 4.556303, delta = 2, p = 1))


overview(mafart_weibull_SC)
plotfit(mafart_weibull_SC, smooth = TRUE)
BIC(mafart_weibull_SC)
mafartparameters_SC <- summary(mafart_weibull_SC)$coefficients
print(mafartparameters_SC)

plot(fitted(mafart_weibull_SC), residuals(mafart_weibull_SC))
abline(h = 0, lty = 2, col = "red")

#mafart ORGANIC
mafart_weibull_SO <- nls(Log_Average_Conc ~ log10_N0 - (Day_numeric/delta)^p, 
                         data =Salm_A_C_1_O, 
                         start = list(log10_N0 = 4.556303, delta = 2, p = 1))


overview(mafart_weibull_SO)
plotfit(mafart_weibull_SO, smooth = TRUE)
BIC(mafart_weibull_SO)
mafartparameters_SO <- summary(mafart_weibull_SO)$coefficients
print(mafartparameters_SO)

plot(fitted(mafart_weibull_SO), residuals(mafart_weibull_SO))
abline(h = 0, lty = 2, col = "red")

#Use lmer to make mixed model to test impact of pond type on log count data using unaveraged aliquots

Salm_lmer <- lmerTest::lmer(log10(Average_Conc) ~ Day * Pond_Type+(1|Pond_Type_Number), Salm_A_C_1)
summary(Salm_lmer)
anova(Salm_lmer)

#add ',type="response" to this emmeans code so it recognizes
#the values are log10 transformed and the responses emmeans provides are exponentiated
#so you dont have to exponentiate by hand
#But you have to calculate % change on your own

emmeans(Salm_lmer, consec ~ Day| Pond_Type,type="response") 
#To calculate % change on your own, Subtract 1 from the factor value (the ratio column in the 
#emmeans output), then multiply by 100

#Conventional Day1/Day0
(0.5166-1)*100 #-48.34 change; decr by ~48%
#Conventional Day2/Day1
(0.3970-1)*100 #-60.3% change; decr by ~60%
#Conventional Day5/Day2
(0.0369-1)*100 #-96.31% change; decr by ~96%


#Organic Day1/Day0
(1.0134-1)*100 #1.34% change; incr by ~1%
#Organic Day2/Day1
(0.6301-1)*100 #-36.99% change; decr by ~37%
#Organic Day5/Day2
(0.0649-1)*100 #-93.51% change; decr by ~93%


#calculate averages

Salm_A_C_bypond <- Salm_A_C_1%>%
  group_by(Day,Pond_Type) %>% 
  mutate(Average_Log_Conc_byPond=log10(mean(Average_Conc)))%>%
  mutate(log_SD_Average_Conc_byPond=sd(Log_Average_Conc)) %>%
  mutate(log_SER_Average_Conc_byPond=(sd(Log_Average_Conc))/(sqrt(3))) %>%
  ungroup


Salm_A_C_bypond <- Salm_A_C_bypond %>%distinct(Day,Pond_Type, .keep_all = TRUE)

Salm_A_C_bypond_fortable <-Salm_A_C_bypond%>%select(Sample_ID, Bacteria,Average_Log_Conc_byPond,log_SD_Average_Conc_byPond)%>%
  mutate(Average_Log_Conc_byPond = round(Average_Log_Conc_byPond, digits = 2))%>%
  mutate(log_SD_Average_Conc_byPond = round(log_SD_Average_Conc_byPond, digits = 2))


#ECOLI_________________________________________________________________

#import E.coli data with select dilutions
Ecoli_select<-read.csv('/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Ecoli_Data_by_Day_forR_07082025_v2.csv')
Ecoli_select$Pond_Type_Number <- paste(Ecoli_select$Pond_Type,Ecoli_select$Pond_Number)

#Make sample ID column
Ecoli_select$Sample_ID_Aliquot<-paste(Ecoli_select$Day,Ecoli_select$Pond_Type_Number,Ecoli_select$Aliquot,sep = '_')
Ecoli_select$Sample_ID<-paste(Ecoli_select$Day,Ecoli_select$Pond_Type_Number,sep = '_')

#calculate average concentration and log average concentration
Ecoli_select1 <- Ecoli_select %>% group_by(Day,Pond_Type_Number,Aliquot,Dilution,Date) %>% mutate(Average_Conc=mean(Final_Conc))%>% mutate(Log_Average_Conc=log10(Average_Conc))%>% ungroup

#get rid of duplicates 
Ecoli_select1<-Ecoli_select1 %>% distinct(Day,Pond_Type_Number,Aliquot,Dilution,Date, .keep_all = TRUE)

#convert sample to character from integer
Ecoli_select1$Pond_Type<-as.character(Ecoli_select1$Pond_Type)
str(Ecoli_select1)

#remove NC samples from Ecoli_select1 to use in mixed model
Ecoli_select1_noNC<-Ecoli_select1%>%  filter(Aliquot!="NC")

##Use lmer to make mixed model to test impact of pond type on log count data using unaveraged aliquots
Ecoli_lmer <- lmerTest::lmer(log10(Average_Conc) ~ Day * Pond_Type+(1|Pond_Type_Number), Ecoli_select1_noNC)
summary(Ecoli_lmer)
anova(Ecoli_lmer)

#add ',type="response" to this emmeans code so it recognizes
#the values are log10 transformed and the responses emmeans provides are exponentiated
#so you dont have to exponentiate by hand
#But you have to calculate % change on your own

emmeans(Ecoli_lmer, consec ~ Day| Pond_Type,type="response") 
#Conventional Day1/Day0
(0.631-1)*100 #-36.9 change; decr by ~37%
#Conventional Day2/Day1
(0.508-1)*100 #-49.2% change; decr by ~49%
#Conventional Day5/Day2
(0.153-1)*100 #-84.7% change; decr by ~85%


#Organic Day1/Day0
(0.870-1)*100 #-13% change; decr by ~13%
#Organic Day2/Day1
(0.843-1)*100 #-15.7% change; decr by ~16%
#Organic Day5/Day2
( 0.149-1)*100 #-85.1% change; decr by ~85%
#remove NC samples, adjust the Sample ID so it has NC attached and then merge with Ecoli_select1_noNC df
Ecoli_select1_NC <- Ecoli_select1%>%
  filter(Aliquot=="NC")
Ecoli_select1_NC$Sample_ID<-paste(Ecoli_select1_NC$Day,Ecoli_select1_NC$Pond_Type_Number,Ecoli_select1_NC$Aliquot,sep = '_')

#Combine data frames
Ecoli_select1v2<-merge(Ecoli_select1_noNC,Ecoli_select1_NC, all=TRUE)

#Arrange so its sorted properly each time
Ecoli_select1v2<-Ecoli_select1v2%>%arrange(Sample_ID_Aliquot)


mycolsEcoli <- rep(c("#E69F00", 
                     "#D55E00", 
                     "#56B4E9",
                     "#0072B2"), 
              length.out = length(unique(Ecoli_select1_noNC$Sample_ID)))

Ecoli_select1_noNC$Day <- factor(Ecoli_select1_noNC$Day, levels = c("DAY0", "DAY1", "DAY2", "DAY5"), 
                         labels = c("Day 0", "Day 1", "Day 2","Day 5"))

Ecoli_select1_noNC<-Ecoli_select1_noNC%>%mutate(Day_numeric=Day)%>%mutate(Day_numeric = str_remove_all(Day, "[^0-9.-]")) %>%
  mutate(Day_numeric = as.numeric(Day_numeric))

LOD.expression <- expression("The Limit of Detection is 1.3"~log [10] ~CFU/mL)

recoli <- ggplot(Ecoli_select1_noNC, aes(x = Day_numeric, y = Log_Average_Conc, 
                                         shape = factor(Pond_Type_Number), 
                                         color = factor(Pond_Type_Number))) + 
  geom_point(size = 18, alpha = 0.9) + 
  labs(title = expression("Enterohemorrhagic" ~ italic("E. coli")),
       x = "Day", 
       y = y.expression) +  
  theme_bw() +
  # Use show.legend = FALSE to keep the dashed line out of your point legend
  geom_hline(yintercept = 1.3, linetype = "dashed", show.legend = FALSE) +
  ylim(1, 5) +
  # Ensure the manual color scale name matches the shape scale (both are NULL/blank)
  scale_color_manual(values = mycolsEcoli) +
  # Optional: force shape values if you don't like the defaults
  scale_shape_manual(values = c(16, 17, 15, 18)) + 
  theme(text = element_text(family = "Times New Roman"),
        strip.text = element_text(size = 40, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank(), # This helps merging
        legend.text = element_text(size = 30),
        axis.text.x = element_text(angle = 0, vjust = 0.3, hjust = 1, size = 40, margin = margin(0,0,20,0)),
        axis.text.y = element_text(angle = 0, vjust = 0.3, hjust = 0.5, size = 40, margin = margin(0,0,0,20)),
        axis.title.y = element_text(size = 50),
        axis.title.x = element_text(size = 50),
        plot.title = element_text(size = 55, face = "italic", hjust = 0.5))

recoli


ggsave(
  filename = "/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Ecoli_contxaxis_poster.tiff", 
  plot = last_plot(),       
  device = "tiff",          
  units = "in", 
  width = 25, 
  height = 15, 
  dpi = 700, 
  compression = "lzw",
  type = "cairo")

#COMPARING LINEAR AND MAFART WEIBULL MODELS
#make a continuous day column

Ecoli_select1_noNC<-Ecoli_select1_noNC%>%mutate(Day_numeric=Day)%>%mutate(Day_numeric = str_remove_all(Day, "[^0-9.-]")) %>%
  mutate(Day_numeric = as.numeric(Day_numeric))

# Fit a model
Ecolimodel <- lm(Log_Average_Conc ~ Day_numeric, data = Ecoli_select1_noNC)

# Generate 4 diagnostic plots (Residuals vs Fitted is first)
plot(Ecolimodel) 
BIC(Ecolimodel)
# To see only the Residuals vs Fitted plot:
plot(Ecolimodel, which = 1)
ggsave("/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Check_linearity/Ecoliresiduals.tiff", units="in", width=20, height=15, dpi=500, compression = 'lzw')

# qqplotplot:
plot(Ecolimodel, which = 2)
ggsave("/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Check_linearity/Ecoliqq.tiff", units="in", width=20, height=15, dpi=500, compression = 'lzw')

#divide Ecoli_noNC in to organic and conventional
Ecoli_select1_noNC_O<-Ecoli_select1_noNC%>%filter(Pond_Type=="Organic")
Ecoli_select1_noNC_C<-Ecoli_select1_noNC%>%filter(Pond_Type=="Conventional")

#linear model
Ecoli_select1_noNC_C_lm <- lm(Log_Average_Conc ~ Day_numeric, data = Ecoli_select1_noNC_C)
BIC(Ecoli_select1_noNC_C_lm)
Ecoli_select1_noNC_O_lm <- lm(Log_Average_Conc ~ Day_numeric, data = Ecoli_select1_noNC_O)
BIC(Ecoli_select1_noNC_O_lm)

#mafart conventional
mafart_weibull_Ecoli_C <- nls(Log_Average_Conc ~ log10_N0 - (Day_numeric/delta)^p, 
                   data =Ecoli_select1_noNC_C, 
                   start = list(log10_N0 = 4.612784, delta = 2, p = 1))


overview(mafart_weibull_Ecoli_C )
plotfit(mafart_weibull_Ecoli_C , smooth = TRUE)
BIC(mafart_weibull_Ecoli_C )
mafartparameters_Ecoli_C  <- summary(mafart_weibull_Ecoli_C )$coefficients
print(mafartparameters_Ecoli_C )
plot(fitted(mafart_weibull_Ecoli_C ), residuals(mafart_weibull_Ecoli_C ))
abline(h = 0, lty = 2, col = "red")

#mafart organic
mafart_weibull_Ecoli_O <- nls(Log_Average_Conc ~ log10_N0 - (Day_numeric/delta)^p, 
                              data =Ecoli_select1_noNC_O, 
                              start = list(log10_N0 = 4.612784, delta = 2, p = 1))


overview(mafart_weibull_Ecoli_O)
plotfit(mafart_weibull_Ecoli_O, smooth = TRUE)
BIC(mafart_weibull_Ecoli_O)
mafartparameters_Ecoli_O  <- summary(mafart_weibull_Ecoli_O)$coefficients
print(mafartparameters_Ecoli_O )
plot(fitted(mafart_weibull_Ecoli_O), residuals(mafart_weibull_Ecoli_O))
abline(h = 0, lty = 2, col = "red")

#get overall averages by organic/conventional status
Ecoli_bypond <- Ecoli_select1 %>% 
  filter(Aliquot!="NC")%>% 
  group_by(Day,Pond_Type) %>% 
  mutate(Average_Log_Conc_byPond=log10(mean(Average_Conc)))%>%
  mutate(log_SD_Average_Conc_byPond=sd(Log_Average_Conc)) %>%
  ungroup

Ecoli_bypond <- Ecoli_bypond %>%distinct(Day,Pond_Type, .keep_all = TRUE)

Ecoli_bypond_fortable <-Ecoli_bypond%>%select(Sample_ID, Bacteria,Average_Log_Conc_byPond,log_SD_Average_Conc_byPond)%>%
  mutate(Average_Log_Conc_byPond = round(Average_Log_Conc_byPond, digits = 2))%>%
  mutate(log_SD_Average_Conc_byPond = round(log_SD_Average_Conc_byPond, digits = 2))


#find average of all EHEC NC samples
Ecoli_select1_NC_overall_avg<-Ecoli_select1_NC%>%mutate(NC_overall_avg=log10(mean(Average_Conc)))


#combine all averaged quantities by pond
Ecoli_bypond_fortable 
Salm_A_C_bypond_fortable 
Lmono_A_C_bypond_fortable 
stacked_df <- bind_rows(Lmono_A_C_bypond_fortable, Salm_A_C_bypond_fortable )

stacked_df <- bind_rows(stacked_df, Ecoli_bypond_fortable )
write.csv(stacked_df, "C:/Users/tmwal/OneDrive/Desktop/stacked_df_TableS2.csv", row.names = TRUE)
