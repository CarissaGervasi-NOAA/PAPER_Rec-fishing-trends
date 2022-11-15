# NOAA Red Snapper Participatory Modeling Project 2022 ##############
# Carissa Gervasi, Mandy Karnauskas, Matt McPherson


library(dplyr)

##
# Load MRIP effort data from online query tool #####################################
##

MRIP = read.csv("../../NOAA Data/PUBLIC/PAPER_Rec-fishing-trends/Raw/MRIP query_effort_81_21_GOM by state_by mode_by area.csv")
summary(MRIP)

MRIP$State = as.factor(MRIP$State)
MRIP$Fishing.Mode = as.factor(MRIP$Fishing.Mode)
MRIP$Fishing.Area = as.factor(MRIP$Fishing.Area)


# The MRIP data is effort in angler trips for each GOM state by fishing area and by mode, 1981-2021


# What are we really interested in? We want to know why recreational anglers fish. Charter and private anglers likely have different motivations. So let's just look at private anglers. Since we are mainly concerned with red snapper, those fish are not being caught from shore. So we really just want the private/rental boat fishing mode. For fishing area, we probably want to remove inland. 

MRIP2 = MRIP %>% 
  filter(Fishing.Mode == "PRIVATE/RENTAL BOAT") %>% 
  filter(Fishing.Area != "INLAND")
MRIP2 = droplevels(MRIP2)
summary(MRIP2)

MRIP2$Fishing.Area2 = as.factor(MRIP2$Fishing.Area)
levels(MRIP2$Fishing.Area2) = list("State" = c("OCEAN (<= 10 MI)", "OCEAN (<= 3 MI)"),
                                   "Federal" = c("OCEAN (> 10 MI)", "OCEAN (> 3 MI)"))

#Remove data points with PSE > 50
library(dplyr)
MRIP2 = MRIP2 %>% 
  filter(PSE < 50)


# Plot everything

MRIP2$effort2 = MRIP2$Angler.Trips/1000000
library(ggplot2)

pdf("Plots/MRIP rec effort by area and state private boats only.pdf", width = 8, height = 5)
ggplot(MRIP2, aes(x=Year, y=effort2, colour=Fishing.Area2)) +
  geom_point() +
  geom_line() +
  facet_wrap(.~State, scales="free_y") +
  theme_bw() +
  labs(x="Year", y="Recreational effort (million angler trips)") +
  scale_x_continuous(breaks=seq(1980,2021,5)) +
  #scale_y_continuous(breaks=seq(0,20,1)) +
  theme(axis.title.x = element_text(vjust=-2),
        axis.title.y = element_text(vjust=5)) + 
  theme(plot.margin = unit(c(0.5,0.8,1,1), "cm"))
dev.off()




# Plot state and federal waters together

MRIP3 = MRIP2 %>% 
  group_by(Year, State) %>% 
  summarize(tottrips = sum(effort2))


pdf("Plots/MRIP rec effort by state private boats only.pdf", width = 8, height = 5)
ggplot(MRIP3, aes(x=Year, y=tottrips)) +
  geom_point() +
  geom_line() +
  geom_smooth() +
  facet_wrap(.~State, scales="free_y") +
  theme_bw() +
  labs(x="Year", y="Recreational effort (million angler trips)") +
  scale_x_continuous(breaks=seq(1980,2021,5)) +
  #scale_y_continuous(breaks=seq(0,20,1)) +
  theme(axis.title.x = element_text(vjust=-2),
        axis.title.y = element_text(vjust=5)) + 
  theme(plot.margin = unit(c(0.5,0.8,1,1), "cm"))
dev.off()

write.csv(MRIP3, "../../NOAA Data/PUBLIC/PAPER_Rec-fishing-trends/Working/Total vessel trips private rec.csv")




# Now let's look at effort just in Florida by wave (private ocean trips only)

FLwave = read.csv("../../NOAA Data/PUBLIC/PAPER_Rec-fishing-trends/Raw/MRIP query_effort_81_21_WFL_private_ocean_by wave.csv")
summary(FLwave)

FLwave$Wave = as.factor(FLwave$Wave)


#Remove data points with PSE > 50
library(dplyr)
FLwave2 = FLwave %>% 
  filter(PSE < 50)


# Plot everything

FLwave2$effort2 = FLwave2$Angler.Trips/1000000
library(ggplot2)

pdf("Plots/MRIP rec effort WFL by wave.pdf", width = 10, height = 5)
ggplot(FLwave2, aes(x=Year, y=effort2)) +
  geom_point() +
  geom_line() +
  geom_smooth() +
  facet_wrap(.~Wave, scales="free_y") +
  theme_bw() +
  labs(x="Year", y="Recreational effort (million angler trips)") +
  scale_x_continuous(breaks=seq(1980,2021,5)) +
  #scale_y_continuous(breaks=seq(0,20,1)) +
  theme(axis.title.x = element_text(vjust=-2),
        axis.title.y = element_text(vjust=5)) + 
  theme(plot.margin = unit(c(0.5,0.8,1,1), "cm"))
dev.off()


pdf("Plots/MRIP rec effort WFL by wave just smoothers.pdf", width = 10, height = 5)
ggplot(FLwave2, aes(x=Year, y=effort2, colour=Wave)) +
  #geom_point() +
  #geom_line() +
  geom_smooth() +
  #facet_wrap(.~Wave, scales="free_y") +
  theme_bw() +
  labs(x="Year", y="Recreational effort (million angler trips)") +
  scale_x_continuous(breaks=seq(1980,2021,2)) +
  #scale_y_continuous(breaks=seq(0,20,1)) +
  theme(axis.title.x = element_text(vjust=-2),
        axis.title.y = element_text(vjust=5)) + 
  theme(plot.margin = unit(c(0.5,0.8,1,1), "cm"))
dev.off()





























# Now remove inshore and see what the trends look like by state


MRIP.off = subset(MRIP, MRIP$Area2 != "Inshore")
MRIP.off = droplevels(MRIP.off)
summary(MRIP.off)

pdf("Plots/MRIP rec effort by area and state no inshore.pdf", width = 8, height = 5)
ggplot(MRIP.off, aes(x=Year, y=effort2, colour=Area2)) +
  geom_point() +
  geom_line() +
  facet_wrap(.~State, scales="free_y") +
  theme_bw() +
  labs(x="Year", y="Recreational effort (million angler trips)") +
  scale_x_continuous(breaks=seq(1980,2022,5)) +
  #scale_y_continuous(breaks=seq(0,20,1)) +
  theme(axis.title.x = element_text(vjust=-2),
        axis.title.y = element_text(vjust=5)) + 
  theme(plot.margin = unit(c(0.5,0.8,1,1), "cm"))
dev.off()



# Let's look at the overall trend without inshore

MRIP.off = MRIP %>% 
  filter(Area2 != "Inshore") %>% 
  group_by(Year) %>% 
  summarize(effort = sum(Angler.Trips))

head(MRIP.off)
str(MRIP.off)
MRIP.off$effort2 = MRIP.off$effort/1000000

library(ggplot2)

pdf("Plots/MRIP rec effort all GOM no inshore.pdf", width = 8, height = 5)
ggplot(MRIP.off, aes(x=Year, y=effort2)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(x="Year", y="Recreational effort (million angler trips)") +
  scale_x_continuous(breaks=seq(1980,2022,2)) +
  scale_y_continuous(breaks=seq(10,40,1)) +
  theme(axis.title.x = element_text(vjust=-2),
        axis.title.y = element_text(vjust=5)) + 
  theme(plot.margin = unit(c(0.5,0.8,1,1), "cm"))
dev.off()



# Next steps

# look at charter trips over time
# compare to license sales
# compare with AL state data effort
# try to get TX state data effort
# compare with red snapper season length
# compare with average gas prices

# Perhaps people are taking more charter trips and less private boats, or rec fishing has gone down with red snapper regulations, or regulations for other species, or gas prices have gotten too high for people






