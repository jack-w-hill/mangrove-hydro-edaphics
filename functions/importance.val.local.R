importance.val.local <- function (z = dataframe)
  
{
  DNAME <- paste(deparse(substitute(z)))
  
  z[z==""] <- NA							# replace empty cells with NA
  
  if ( any(is.na(z)) )
    stop("Empty cells are not allowed in the data frame ", DNAME)
  
  if ( dim(z)[2] != 5 )
    stop(DNAME, " must have 5 columns in this order: Point, Quarter, Species, Distance, DBH")
  
  ## rename the columns
  colnames(z) <- c('Point', 'Quarter', 'Species', 'Distance', 'DBH')
  
  ## there should be no vacant quarters, but check again:
  if ( any(is.na(z[['Distance']])) )
    stop(DNAME, "Empty cells are not allowed. Recheck data.")
  
  
  n = nlevels(factor(z[['Point']]))		# number of sample points along transect
  quarters = length(factor(z[['Point']]))	# number of observations
  
  if ( quarters != 4*n | max(table(z$Point)) != min(table(z$Point)) )
    stop("Some sample point does not have exactly four quarters in data frame ", DNAME)
  
  if ( length(table(z$Quarter)) != 4 | max(table(z$Quarter)) != min(table(z$Quarter)) ) 
    warning("All sample points in data frame ", DNAME, " do not use the same four quarter names.")
  
  #species = factor(z[['Species']])			# list species observed
  
  # TS MODIFICATION #
  # The commented-out code only allowed for calculations for species observed in
  # the transect, not absent species that were observed in other transects.
  # Re-writing so it is easier to iterate over multiple transects with different
  # assemblages
  species = z[['Species']]			# list species observed
  
  ## do the analysis
  BA = pi * z[['DBH']]^2 /4					# basal area for each tree = pi * d^2 / 4
  
  r = mean(z[['Distance']])					# mean distance
  absDensity = 1 / r^2						# absolute density per m^2
  absDensityPerHA = 10000 * absDensity		# absolute density per ha
  
  sumBA = tapply(BA,species,sum)		# total basal area by species
  sumBA[is.na(sumBA)] = 0
  
  meanBA = tapply(BA,list(species),mean)		# mean basal area by species
  meanBA[is.na(meanBA)] = 0
  
  prop = table(z[['Species']]) / quarters		# proportion of each species in all observations
  relDensity = round(100 * prop, 2)			# relative density by species (as %)
  
  treesPerHA = prop * absDensityPerHA			# number of trees per species per ha
  absCover = meanBA * treesPerHA/10000		# cover (in m^2) by species per ha
  totalCover = sum(absCover)					# total cover (in m^2) by all species per ha
  
  relCover = round(100 * absCover/totalCover, 2)	# relative cover as a percent
  
  ## create matrix of points by species: 1 = species is present at point, 0 = otherwise
  ptXsp = ceiling(table(z[['Point']], z[['Species']]) / 4) 
  
  ## number of points (not quarters) at which each species is observed
  absObservations = colSums(ptXsp, na.rm = FALSE, dims = 1) 
  
  relFrequency = round(100 * absObservations / sum(absObservations), 2)	# relative frequency by species
  
  importance = relDensity + relCover + relFrequency		# abs importance by species
  relImportance = 100 * importance/sum(importance)		# rel importance by species
  
  absDensity <- absDensityPerHA * table(z$Species) / quarters
  
  ## bind the results into a data frame
  results <- as.data.frame(rbind(relDensity, relCover, relFrequency, importance, relImportance, absDensity))
  
  ## transpose the data frame so that each row represents a species
  results <- as.data.frame(t(results))
  
  ## make the column names
  colnames(results) <- c(" Rel Density", " R Cover", " R Freq", " Importance", " R Import", " Abs Density")
  
  ## print results in decreasing order of importance
  
  cat("Number of sample points: n =", n, "\n")    
  cat("Overall Absolute Density per Hectare (Cottam & Curtis):", format(round(absDensityPerHA, 2), nsmall = 2), "\n", "\n")    
  
  # TS modification - remove formatting and rounding of output table
  return(results[order(-importance), ])	
}
