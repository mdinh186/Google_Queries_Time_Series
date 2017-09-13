library(data.table)

ultimatePai = function (test, uniqueBool = T){
  #modify incident name of testprob, test

  test$Incident.Name = gsub(x = test$Incident.Name, pattern = " ", replacement = ".")
  if (uniqueBool == T) {
    test = unique(test[,c("Week","Incident.Name","Weekday","Hour","Month","Day", "Year","Black","Vacant.Housing","Mean.Income","Income.0.10000","Leave.to.Work.12.5am", "Weekday Density", "Weekly Density","Spatial Density", "Hourly 2d Density", "Hourly 1d Density", "Weekday1dDensity", "Week1dDensity","WeekDens","OverallWeekDens", "Cell", "Cost")])
    testProb = readRDS("testProb2.rds")
    testProb$Incident.Name = gsub(x = testProb$Incident.Name, pattern = " ", replacement = ".")
    crimetype = colnames(testProb)[18:35]
    crimecounts = readRDS('predictionsListUniqueGBM.rds')
  }else{
    testProb = readRDS("testProb.rds")
    testProb$Incident.Name = gsub(x = testProb$Incident.Name, pattern = " ", replacement = ".")
    crimetype = colnames(testProb)[24:41]
    crimecounts = readRDS('predictionsListNonUnique.rds')
    }
  summary = data.table(Cell = unique(test$Cell))
  paitable = data.table(matrix(0, ncol=7))
  for (i in 1:18){
    if (uniqueBool == T){
      dat = testProb[Incident.Name == crimetype[i],][order(Cell),17:35]
    } else{
      dat = testProb[Incident.Name == crimetype[i],][order(Cell),c(22,24:41)]
      }

    dat = data.table(Cell = dat$Cell, Prob = dat[[i+1]])
    crimeCount = data.table(crimecounts[[i]])# with 7 columns
    MaxCount = data.table(crimeCount[,-1] * dat$Prob)  # change column to -1 when have new data
    Cost = test[Incident.Name == crimetype[i], Cost]
    currCrimeCost = MaxCount * Cost
    currCrimeCost$Cell = dat$Cell
    name = head(colnames(currCrimeCost),-1)
    colnames(currCrimeCost) = c(paste0(crimetype[i], ".",name), "Cell")
    setcolorder(currCrimeCost,c("Cell", names(currCrimeCost)[1:length(currCrimeCost)-1]))
    currCrimeCost2 = data.table(aggregate(currCrimeCost[,-1], by = list(Cell = currCrimeCost$Cell), FUN = sum))
    summary = merge(summary, currCrimeCost2, all.x = TRUE, by= "Cell")
    summary[is.na(summary)] =0
  }
  colnames(paitable) = name
  #caculate pai for each model:
    for (i in 2:8){
      currModel = data.table(Cell = summary$Cell)
      colidx = seq(i, ncol(summary), by =7)
      for (idx in colidx){
        currModel= cbind(currModel, summary[[idx]])
      }
      colnames(currModel) = c("Cell",crimetype)
      currModel$TotalCost= apply(currModel[,c(2:19)], 1, sum)
      currModel = currModel[order(-TotalCost),]
      nbflagged = ceiling(nrow(currModel) * 0.01)
      flagged_cost = sum(currModel[1:nbflagged, TotalCost])
      total_cost = sum(currModel[,TotalCost])
      fac1 = flagged_cost/total_cost
      fac2 = nrow(currModel)/nbflagged
      pai = fac1 *fac2
      paitable[[i-1]] = pai
    }
  return (paitable)
}


#test = unique(test[,c("Week", "Incident.Name", "Hour", "Month", "Day", "Year", "Weekday", "WeekDens", "OverallWeekDens","Weekly Density","Spatial Density","Hourly 2d Density","Hourly 1d Density", "Weekday Density", "Weekday1dDensity", "Week1dDensity", "Cell", "Cost")])
#probU = readRDS("testProb.rds")

test = readRDS("test.rds")
train = readRDS('train.rds')


ultimatePai(test = test, uniqueBool = T)
