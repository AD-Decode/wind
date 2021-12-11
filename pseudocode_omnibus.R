
#https://cran.r-project.org/web/packages/afex/vignettes/afex_anova_example.html 

sink(paste(outpath,'Probe_d5_tukeytable.txt', sep='')) 



lm_Probe_d8 <- lm(Probe_d8 ~ Genotype*Sex, geno_combined_FA) 

summary(lm_Probe_d5)  #gives model test p value 

anova(lm_Probe_d5) #gives significance of factors 

summary(glht(lm_Probe_d5, emm(pairwise ~ Genotype*Sex, contr="sidak"))) #consider save to file #a whole slew of pairwise comparisons 

summary(glht(lm_Probe_d5, emm(pairwise ~ Genotype*Sex, contr="sidak", by=NULL)))  

#emmeans(lm_Probe_d5, ~ Genotype*Sex, contr="tukey") 



summary(as.glht(pairs(Probe_d5.emm), by = NULL)) 

summary(as.glht(pairs(Probe_d5.emm)), test=adjusted("free")) 



mypairs<-pairs(Probe_d5.emm) #_consider save to file 

pairs(Probe_d5.emm, by="Genotype") 

pairs(Probe_d5.emm, by="Sex") 

sink() 



#omnibus testing 

#https://cran.r-project.org/web/packages/emmeans/vignettes/confidence-intervals.html#adjust 

sink() 

sink(paste(outpath, "Winding_ProbeD5_omnibus.txt")) 

test(mypairs, joint = TRUE, adjust = "sidak") 

joint_tests(ref_grid(lm_Probe_d5)) 

joint_tests(ref_grid(lm_Probe_d5), by = "Sex", adjust = "sidak") 

joint_tests(ref_grid(lm_Probe_d5), by = "Genotype",adjust = "sidak" ) 

sink()