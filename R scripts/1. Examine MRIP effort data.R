
# NOAA Red Snapper Participatory Modeling Project 2022 ##############
# Carissa Gervasi, Mandy Karnauskas, Matt McPherson


##
# Load MRIP effort data from online query tool #####################################
##

MRIP = read.csv("Data/Raw//MRIP query_effort_81_21_GOM by state_private_by area.csv")
summary(MRIP)

# The MRIP data is effort in angler trips for each GOM state by fishing area, private/rental boats only, 1981-2021

# First let's look at the overall trend

library(dplyr)

MRIP.all = MRIP %>% 
  group_by(Year) %>% 
  summarize(effort = sum(Angler.Trips))

head(MRIP.all)
str(MRIP.all)
MRIP.all$effort2 = MRIP.all$effort/1000000

library(ggplot2)

pdf("Plots/MRIP rec effort all GOM.pdf", width = 8, height = 5)
ggplot(MRIP.all, aes(x=Year, y=effort2)) +
  geom_point() +
  geom_line() +
  theme_classic() +
  labs(x="Year", y="Recreational effort (million angler trips)") +
  scale_x_continuous(breaks=seq(1980,2022,2)) +
  scale_y_continuous(breaks=seq(10,40,1)) +
  theme(axis.title.x = element_text(vjust=-2),
        axis.title.y = element_text(vjust=5)) + 
  theme(plot.margin = unit(c(0.5,0.8,1,1), "cm"))
dev.off()



# Now let's break down the data by state and by area

MRIP$Area2 = as.factor(MRIP$Fishing.Area)
levels(MRIP$Area2)
levels(MRIP$Area2) = list("Inshore" = c("INLAND"), "State" = c("OCEAN (<= 10 MI)","OCEAN (<= 3 MI)"), "Federal" = c("OCEAN (> 10 MI)","OCEAN (> 3 MI)"))


MRIP$State = as.factor(MRIP$State)
levels(MRIP$State)

MRIP$effort2 = MRIP$Angler.Trips/1000000


pdf("Plots/MRIP rec effort by area and state.pdf", width = 8, height = 5)
ggplot(MRIP, aes(x=Year, y=effort2, colour=Area2)) +
  geom_point() +
  geom_line() +
  facet_wrap(.~State, scales="free_y") +
  theme_classic() +
  labs(x="Year", y="Recreational effort (million angler trips)") +
  scale_x_continuous(breaks=seq(1980,2022,5)) +
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
  theme_classic() +
  labs(x="Year", y="Recreational effort (million angler trips)") +
  scale_x_continuous(breaks=seq(1980,2022,5)) +
  #scale_y_continuous(breaks=seq(0,20,1)) +
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






