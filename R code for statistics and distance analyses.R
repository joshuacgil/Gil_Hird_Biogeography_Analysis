####R Code used for statisical analyses####

##Packages required for the following code##
library(vegan); packageVersion("vegan")
library(phyloseq); packageVersion("phyloseq")
library(Imap);packageVersion('Imap')
library(sp);packageVersion("sp")
library(gstat);packageVersion("gstat")
library(spaa);packageVersion("spaa")
library(Imap);packageVersion("Imap")
library(ggplot2); packageVersion("ggplot2")



##Alpha Diversity##

#Calcualte Shannon diveristy index for every sample
shannon = estimate_richness($phyloseq_object, split=TRUE, measures = 'shannon')

#add alpha diversity to phyloseq object's metadata file
sample_data($phyloseq_object)$shannon <- estimate_richness($phyloseq_object, measures='shannon')

meta = sample_data($phyloseq_object)[order(row.names(sample_data($phyloseq_object))),]

#Calcaulte differences between ultiple groups
kruskal.test(as.numeric(unlist(meta$shannon)), as.factor(meta$Varaible_of_Interest))

#Wilcox rank-sum calcaulats all pairwise comparisons and generates a matrix showing all compariosns
pairwise.wilcox.test(as.numeric(unlist(meta$shannon)), as.factor(meta$Varaible_of_Interest),  p.adjust.method='fdr')



##Beta Diverisity##

#Generating the PCoA ordination, Options Used are: 'wunifrac', 'unifrac', and 'bray'

#Example, Weighted Unifrac
set.seed(1)
phyloseq_wunifrac<-ordinate($phyloseq_object, "PCoA", "wunifrac")
ordination <-plot_ordination($phyloseq_object, phyloseq_wunifrac, color="Variable_of_Interest", label=NULL, shape=NULL)+ ggtitle("PCoA Ordination") + geom_point(size=5)
ordination

set.seed(1)

#Methods = wunifrac, uunifrac, bray
weighted_nova <- phyloseq::distance($phyloseq_object, method = "wunifrac")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data($phyloseq_object))

# Adonis/PERMANOVA test
adonis(weighted_nova ~ Variable_of_interest, data = sampledf, permutations =999)



##Distance Decay Analaysis##

#Run The functions 'ReplaceLowerOrUpperTriangle' and 'GeoDistanceInMetresMatrix'. These are used to create a distance matrix from longitude and latitude data. 
ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.
  
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.
  
  GeoDistanceInMetres <- function(g1, g2){
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$name > g2$name, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="m")))
    }
    return(mapply(DistM, g1, g2))
  }
  
  n.geopoints <- nrow(df.geopoints)
  
  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$name <- 1:n.geopoints
  
  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("name", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})
  
  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
  
  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name
  
  return(mat.distances)
}



##The following merges the beta diversity distance matrix and the geographic distance matrix used to make distance decay plots 

#calculate beta diveristy distance matrix 
weighted_nova <- phyloseq::distance($phyloseq_object, method = "wunifrac")
unweighted_nova <- phyloseq::distance($phyloseq_object, method = "unifrac")
bray_nova <- phyloseq::distance($phyloseq_object, method = "bray")

#Metadata file has to have the 'Latitude' and 'Longitude' in separate columns
sd = data.frame(sample_data($phyloseq_object))
sd <- sd[,c("Latitude","Longitude")]
sd[,c('Latitude')]
row.names(sd)
name <- row.names(sd)
lat <- sd[,c("Latitude")]
lon <- sd[,c("Longitude")]

df <- data.frame(name, lat, lon)
df

#converts lat and lon to numerics 
df$lat <- as.numeric(as.character(df$lat))
df$lon <- as.numeric(as.character(df$lon))

#generates a distance matrix
dist.mat <- GeoDistanceInMetresMatrix(df)
dist.mat

#convert Distance matrix into dataframe, use: 'weighted_nova', 'unweighted_nova', or 'bray_nova' 
dist_metric_df <-dist2list(weighted_nova) 
dist_metric_df


#Optional, convert the distance matrix into Kilometers
dist.mat <-GeoDistanceInMetresMatrix(df)/1000
dist.mat

#clean up entire distance matrix
dist.mat
geo_dist_object <- as.dist(dist.mat)
geo_dist_object

#converts the total distances into a data frame that can be joint to the other data frame 
geo_dist_df <- dist2list(geo_dist_object)
geo_dist_df

#mergeing the dataframes by appending the column from the beta dist
joined.variogram <- geo_dist_df
joined.variogram$beta_diversity = dist_metric_df$value #in the dist_metric_df bray curtis data is saved as "Value"

# remove any null data rows
variogram.df <- joined.variogram
variogram.df=na.omit(variogram.df)
variogram.df

#remove the zeros (i.e. self comparions) in the distances matrices
variogram.df <- variogram.df[variogram.df$beta_diversity != 0,]
variogram.df

#Make the distance decay plot with generalized additive model regression line. 
ggplot(data = variogram.df, aes(x = value, y = beta_diversity, colour=NULL)) + 
  geom_point(size=1) + 
  facet_grid(NULL) + 
  geom_point(color = "grey75") +
  geom_smooth(method = "gam", formula=y~s(x), show.legend = TRUE)+
  labs(title = "Distance Decay weighted UniFrac diveristy")+
  labs(x="Distance in kilometers")+
  labs(y="Weighted UniFrac Distance Dissimilarity")+
  theme_bw()


##Mantel test to test for signifcant correlation between geography and beta diveristy##
mantel(weighted_nova, geo_dist_object, method="spearman", permutations=999 )
mantel(unweighted_nova, geo_dist_object, method="spearman", permutations=999 )
mantel(bray_nova, geo_dist_object, method="spearman", permutations=999 )

