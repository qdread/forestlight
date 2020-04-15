#from Kitijama, 2005 Annals of Botany, using datathief 1.7 to digitize data
library(tidyverse)
library(broom)
library(lme4)

lai <- read_csv('/Users/jgradym/Google_Drive/ForestLight/data/data_forplotting/LAI_Depth.csv')
pfd <- read_csv('/Users/jgradym/Google_Drive/ForestLight/data/data_forplotting/PFD_LAI.csv')

ggplot(data = lai, 
       aes(x = Depth, y = LAI, fill = Species, color = Species )) +
  geom_point(size = 5, shape = 21, color = "black") + theme_plant  +
  scale_y_log10(labels = signif) +
  scale_x_log10(labels = signif, name = "Crown Depth (m)") +
  geom_smooth(method = "lm", alpha = 0.2, aes(fill = Species, color = Species))

p <- ggplot(data = lai %>% filter(Species != "Cecropia"), 
            aes(x = Depth, y = LAI, fill = Species, color = Species )) +
  geom_point(size = 4.5, shape = 21, color = "black") + theme_plant  +
  scale_y_log10(limits = c(0.6, 10),labels = signif, breaks = c(1, 3, 10)) +
  scale_x_log10(limits = c(0.2, 15), breaks = c( 0.3, 1, 3, 10), labels = signif, name = "Crown Depth (m)") +
  geom_smooth(method = "lm", alpha = 0.2, aes(fill = Species, color = Species))
p


p2 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(p2)
pdf(file.path(gdrive_path,'Figures/Supplementals/LAI/lai_depth.pdf'))
grid.draw(p2)
dev.off()
system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path,'Figures/Supplementals/LAI/lai_depth.pdf'), 
                    file.path(gdrive_path,'Figures/Supplementals/LAI/lai_depth.pdf')) 
)

# with species
lm1 <- lm(log(LAI) ~ log(Depth)+ Species, data = lai)
summary(lm1) # 0.351, r2 = 0.99

# without species
lm2 <- lm(log(LAI) ~ log(Depth)+ Species, data = lai)
summary(lm2) # 0.54, r2 = 0.38

# calculated separately
lm3 <- lm(log(LAI) ~ log(Depth), data = lai2)
summary(lm3)
lm_each <-  lai %>%
  nest(-Species) %>% #group variable
  mutate(
    fit = map(data, ~ lm(log(LAI) ~ log(Depth), data = .x)),
    tidied = map(fit, tidy)
  ) %>%
  unnest(tidied)%>%
  filter(term != '(Intercept)')
lm_each
mean(lm_each$estimate) #0.331


# Mixed Model
lai2 <- lai %>% filter(Species != "Cecropia") #Cecropia is weird

# all 5 species
mlm <- lmer(log(LAI) ~ log(Depth) + (log(Depth) |Species), lai) 
mlm # 0.33 sloep

#no Cecropia 
mlm2 <- lmer(log(LAI) ~ log(Depth) + (log(Depth) |Species), lai2) 
mlm2 #0.37 slope

### PFD analysis
ggplot(data = pfd,aes( x = LAI, y = PFD, fill = Species, color = Species )) +
  geom_point(size = 6, shape = 21, color = "black") + theme_plant  +
  scale_y_log10() +
  scale_x_log10(limits = c(0.3, 10)) +
  geom_smooth(method = "lm", alpha = 0.2, aes(fill = Species, color = Species))

#semi-log
p <- ggplot(data = pfd %>% filter(Species != "Cecropia"),
       aes( x = LAI, y = PFD, fill = Species, color = Species )) +
  geom_point(size = 4.5, shape = 21, color = "black") + theme_plant  +
  scale_y_log10(limits = c(2.5, 110), name = "% PFD Transmission") +
  scale_x_continuous(limits = c(-0.3, 8)) +
  geom_smooth(method = "lm", alpha = 0.2, aes(fill = Species, color = Species))
p
pdf(file.path(gdrive_path,'Figures/Supplementals/LAI/lai_pfd.pdf'))
p
dev.off()

p2 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(7,"cm"))
grid.newpage()
grid.draw(p2)
pdf(file.path(gdrive_path,'Figures/Supplementals/LAI/lai_pfd.pdf'))
grid.draw(p2)
dev.off()
system2(command = "pdfcrop", 
        args    = c(file.path(gdrive_path,'Figures/Supplementals/LAI/lai_pfd.pdf'), 
                    file.path(gdrive_path,'Figures/Supplementals/LAI/lai_pfd.pdf')) 
)

# all 5 species
mlm_pdf <- lmer(log(PFD) ~ LAI + (LAI |Species), pfd) 
mlm_pdf # 0.33 sloep

#no cecropia
pfd2 <- pfd %>% filter(Species != "Cecropia")
mlm_pdf2 <- lmer(log(PFD) ~ LAI + (LAI |Species), pfd2) 
mlm_pdf # 0.33 sloep

#no Cecropia 

mlm2_pdf <- lmer(log(PFD) ~ LAI + (LAI |Species), pfd2) 
mlm2_pdf  #0.37 slope

# pfd all
lm_pfd <- lm(log(PFD) ~ LAI, pfd)
lm_pfd
lm2_pfd <- lm(log(PFD) ~ LAI + Species, pfd)
lm2_pfd

lm3_pfd <- lm(log(PFD) ~ LAI + Species, pfd2)
lm3_pfd

lm_each_pfd <-  pfd %>%
  nest(-Species) %>% #group variable
  mutate(
    fit = map(data, ~ lm(log(PFD) ~ LAI, data = .x)),
    tidied = map(fit, tidy)
  ) %>%
  unnest(tidied)%>%
  filter(term != '(Intercept)')
lm_each_pfd
mean(lm_each_pfd$estimate) #0.331


lm_each_pfd2 <-  pfd2 %>%
  nest(-Species) %>% #group variable
  mutate(
    fit = map(data, ~ lm(log(PFD) ~ LAI, data = .x)),
    tidied = map(fit, tidy)
  ) %>%
  unnest(tidied)%>%
  filter(term != '(Intercept)')
lm_each_pfd2
mean(lm_each_pfd2$estimate) #0.331
