rm(list=ls()) #clean memory
gc()          #collect garbage

############################################# --- begin user data --- #############################################
datadir="D:/data/om_panels"   #where is the data located
codedir="B:/GitHub/Platform-Drosophila/Platform-Optogenetics" #location of other R or Rmd files used in this script, normally location of this script
htmlname="test.html" #filename for HTML evaluation sheet
groupfilename="test.txt" #filename for text file with datafiles assigned to experimental groups
############################################## --- end user data --- ##############################################

library(ggplot2)

#################################################################### Functions ####################################################


### Function to downsample the data using the approx() function (for data with period/time jitter)
downsampleapprox <- function(data, experimentDuration, NofPeriods, NofDatapoints) {
  
  # create vectors for ball position, time and period number
  sensRotYDownsampled <- vector(mode = "numeric")
  timeDownsampled = seq(0, (as.numeric(as.character(experimentDuration))*1000)-50, 50) #create a time vector (may not be necessary)
  periodDownsampled <- rep(1:NofPeriods, each = 400) #create period numbers

  # downsample ball position
  for (index in 1:NofPeriods){
    p=round(approx(subset(data$SensRotY*10000, data$CurrentTrial==index), n=table(periodDownsampled)[index])$y)/10000
    p[p < 0]=p[p < 0] + 360 #in case this has generated negative values
    p[p > 360]=p[p > 360] - 360 #in case this has generated too high values
    sensRotYDownsampled=c(sensRotYDownsampled, p) #Concatenate the traces from each period
  }
  
  # bind the downsampled vectors into one dataframe
  dataDown <- data.frame("time" = timeDownsampled, "period" = periodDownsampled, "a_pos" = sensRotYDownsampled)
  
  # return the downsampled data
  return(dataDown)
}




#initialize some values
NofPeriods = 48 #3 times 8 patterns l+r
experimentDuration = 960 #3 times 8 patterns l+r for 20s each
NofDatapoints = 19200 #3 times 8 patterns l+r for 20s each at 20Hz
rot_amount <- vector(mode = "numeric")
sum_rot <- vector(mode = "numeric")


#Load the file that contains all the saved files
tested_flies <- read.table(paste(datadir,"/",groupfilename, sep = ""), quote="\"", comment.char="#")$V1 #load file with fly names
NofFlies=length(tested_flies)
setwd(datadir)
grp_rot <- data.frame(matrix(ncol = NofFlies, nrow = 16)) #a dataframe to put all rotatory data of each fly

for (f in 1:NofFlies) {
  rawdata <- read.csv(tested_flies[f]) #read the raw data from the fly
  
  #Generate approx. period numbers
  tempPeriods<-rep(1:NofPeriods, each = trunc(nrow(rawdata)/48))
  rawdata$CurrentTrial[1:length(tempPeriods)] <- tempPeriods
  rawdata$CurrentTrial[(length(tempPeriods) + 1):nrow(rawdata)] <- 48  # fill any remaining values with the last period
  
  lodata <- downsampleapprox(rawdata, experimentDuration, NofPeriods, NofDatapoints) #downsample the data to 20Hz
  rotresp <- diff(lodata$a_pos) #calc the difference between two angle values
  rotresp[rotresp > 180]=rotresp[rotresp > 180]-360   #take the smaller angles
  rotresp[rotresp < -180]=rotresp[rotresp < -180]+360 #take the smaller angles
  lodata$rotresp[1] = 0
  lodata$rotresp[-1] = rotresp #add the rotatory response to the dataframe  

  #sum the rotatory responses by period
  for (period in 1:NofPeriods) {
    rot_amount[period]=sum(lodata$rotresp[lodata$period==period])  
  }
  
  #sum the rerponses in each repeated period
  for (x in 1:16) {
    sum_rot[x]=rot_amount[x]+rot_amount[x+16]+rot_amount[x+32]
  }
  barplot(sum_rot)  
  grp_rot[f] = sum_rot
}  # done with all the flies

grp_rot[c(2,4,6,8,10,12,14,16), ]=-grp_rot[c(2,4,6,8,10,12,14,16), ]
optomotor <- data.frame(matrix(ncol = NofFlies, nrow = 8)) #a dataframe for the averaged responses
for (i in 1:8) {
  row_idx1 <- 2*i - 1
  row_idx2 <- 2*i
  optomotor[i,] <- (grp_rot[row_idx1,] + grp_rot[row_idx2,]) / 2
}

#plot results:
optomotor$wavelength=c(4.8,5.8,7.2,9.7,15,20,24,40) 
optomotor$stripes=c(75,62,50,37,24,18,15,9)
optomotor$mean <- rowMeans(optomotor)
optomotor$se <- apply(optomotor[, 1:NofFlies], 1, function(x) sd(x)/sqrt(length(x)))


# Create a plot with error bars
# Create a new dataframe for plotting with row numbers as x-axis
plot_df <- data.frame(
  row_num = 1:nrow(optomotor),
  mean = optomotor$mean,
  se = optomotor$se,
  labels2 = optomotor$wavelength,
  labels1 = optomotor$stripes
)

# Create the plot
print(ggplot(plot_df, aes(x = labels1, y = mean)) +
  # Add horizontal dashed grey line at y=0 (in the background)
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.7) +# Add smooth curve first (so it appears behind points)
  geom_smooth(method = "loess", span = .75, se = FALSE, color = "blue", size = 1) +        
  #now the points and error bars
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_x_continuous(
      name = "Number of stripes",
      sec.axis = sec_axis(
        transform = ~.,  # Identity transformation - same scale
        name = "Pattern Wavelength",
        breaks = plot_df$labels1,
        labels = plot_df$labels2
      )
    ) +
    labs(
      title = "Rotatory Response (+/- SEM)",
      y = "Cumulative Rotatory Response [deg]"
    ) +
    theme_classic() +
    theme(
      axis.title.x.top = element_text(color = "blue"),
      axis.text.x.top = element_text(color = "blue")
    ))