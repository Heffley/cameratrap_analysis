library(MuMIn)
library(unmarked)
library(dplyr)

###################START HERE IF WANTING TO USE DIRECTLY THE DETECTION MATRIX #############
###read directly the detection matrix RDS to avoid every step of this script up to here###
detection_matrix  <- readRDS(gzcon(url("https://github.com/tgelmi-candusso/cameratrap_analysis/raw/main/detection_matrix_Scarborough.rds")))

#### COVARIATES ########

##GENERATE COVARIATE dataframes for the model , make sure to readapt the site_names AND add human/dog presence####
urlfile500="https://raw.githubusercontent.com/tgelmi-candusso/cameratrap_analysis/main/cov_500.csv"
urlfile1000="https://raw.githubusercontent.com/tgelmi-candusso/cameratrap_analysis/main/cov_1000.csv"
urlfile2000="https://raw.githubusercontent.com/tgelmi-candusso/cameratrap_analysis/main/cov_2000.csv"
urlfilehumans="https://raw.githubusercontent.com/tgelmi-candusso/cameratrap_analysis/main/human_dog_df.csv"

human_dog_df <- read.csv(urlfilehumans) %>% 
  select(-1) %>% 
  select(site_name, total_freq_humans ,total_freq_dogs )

b500 <- read_csv(urlfile500)%>%
  mutate(site_name = gsub("_", "", site_name))%>%
  mutate(site_name = gsub("TUW0", "TUW", site_name))
b500 <- left_join(b500, human_dog_df, by="site_name")%>%
  dplyr::filter(site_name %in% unique(detection_matrix$deer$site_name)) ##filter those for relevant for the analysis

b1000 <- read_csv(urlfile1000)%>%
  mutate(site_name = gsub("_", "", site_name))%>%
  mutate(site_name = gsub("TUW0", "TUW", site_name))
b1000 <- left_join(b1000, human_dog_df, by="site_name")%>%
  dplyr::filter(site_name %in% unique(detection_matrix$deer$site_name)) ##filter those for relevant for the analysis


b2000 <- read_csv(urlfile2000)%>%
  mutate(site_name = gsub("_", "", site_name))%>%
  mutate(site_name = gsub("TUW0", "TUW", site_name))
b2000 <- left_join(b2000, human_dog_df, by="site_name")%>%
  dplyr::filter(site_name %in% unique(detection_matrix$deer$site_name)) ##filter those for relevant for the analysis



##call occupancy covariates
cov<- b1000 %>% select(-1, -BUFF_DIST, -SHAPE_Length, -ORIG_FID, -SHAPE_Area)
cov <- as.data.frame(scale(cov))

##call detection covariate matrix here if using
##det_list <- list(season = det_covs)

# setting up for occupancy for deer
y <- detection_matrix$deer[ , 2:54]

siteCovs_500 <- as.data.frame(b500)
siteCovs_500 <- siteCovs_500[,2:28]

siteCovs_1000 <- as.data.frame(b1000)
siteCovs_1000 <- siteCovs_500[,2:28]

siteCovs_2000 <- as.data.frame(b2000)
siteCovs_2000 <- siteCovs_500[,2:28]

umf_deer_500 <- unmarkedFrameOccu(y = y, siteCovs = siteCovs_500)
umf_deer_1000 <- unmarkedFrameOccu(y = y, siteCovs = siteCovs_1000)
umf_deer_2000 <- unmarkedFrameOccu(y = y, siteCovs = siteCovs_2000)

##SINGLE SPECES
y_list_c <- list(coyote = as.matrix(detection_matrix$coyote %>% select(-1)))
coyote <- unmarkedFrameOccu(y = (detection_matrix$coyote%>%select(-1)),
                            siteCovs = cov)
#obsCovs = det_list
#mdata <- coyote
## in order to change deer model with buffer change swap mdata variable with umf_deer_500, umf_deer_1000, or umf_deer_2000
mdata <- umf_deer_2000


##single covariate comparison

fit_null <- occu(formula = ~ 1
                 ~ 1,
                 data = mdata)


fit_null <- occu(formula = ~1
                      ~1,
                 data = mdata)

fit_LFT <- occu(formula = ~1
                     ~LFT_dist,
                data = mdata)

fit_H2O <- occu(formula = ~1
                     ~H2O_dist,
                data = mdata)

fit_WV <- occu(formula = ~1
                    ~WV_dist,
               data = mdata)

fit_MV <- occu(formula = ~1
                    ~MV_dist,
               data = mdata)

fit_WVF <- occu(formula = ~1
                     ~WVF_dist,
                data = mdata)

fit_WVO <- occu(formula = ~1
                     ~WVO_dist,
                data = mdata)

fit_built <- occu(formula = ~1
                  ~built,
                  data = mdata)

fit_DEM_median <- occu(formula = ~1
                            ~DEM_median,
                            
                            data = mdata)

fit_DEM_mean <- occu(formula = ~1
                          ~DEM_mean,
                          
                          data = mdata)

fit_NDVI_median <- occu(formula = ~1
                             ~NDVI_median,
                             
                             data = mdata)

fit_NDVI_mean <- occu(formula = ~1
                           ~NDVI_mean,
                           
                           data = mdata)

fit_POP_median <- occu(formula = ~1
                            ~POP_median,
                            
                            data = mdata)

fit_POP_mean <- occu(formula = ~1
                          ~POP_mean,
                          
                          data = mdata)

fit_WVO_PA <- occu(formula = ~1
                        ~WVO_PA,
                        
                        data = mdata)

fit_WVF_PA <- occu(formula = ~1
                        ~WVF_PA,
    
                        data = mdata)

fit_MV_PA <- occu(formula = ~1
                       ~MV_PA,
                       
                       data = mdata)

fit_FD_PA <- occu(formula = ~1
                       ~Fdec_PA,
                       
                       data = mdata)

fit_FM_PA <- occu(formula = ~1
                       ~Fmix_PA,
                       
                       data = mdata)

fit_FC_PA <- occu(formula = ~1
                       ~Fcon_PA,
                       
                       data = mdata)
fit_cor <- occu(formula = ~1
                     ~corridor,
                     
                     data = mdata)
fit_hum <- occu(formula = ~1
                     ~total_freq_humans,
                     
                     data = mdata)
fit_dog <- occu(formula = ~1
                     ~total_freq_dogs,
                     
                     data = mdata)

fit <- fitList(fit_null, fit_LFT, fit_H2O, fit_WV, fit_MV, fit_WVF, fit_WVO,
               fit_built, fit_DEM_median, fit_DEM_mean, fit_NDVI_median,
               fit_NDVI_mean, fit_POP_median, fit_POP_mean, fit_WVO_PA, fit_WVF_PA, 
               fit_MV_PA, 
               fit_FC_PA, 
               fit_FM_PA, fit_FD_PA,
               fit_cor, fit_hum, fit_dog)
modSel(fit)

coy_SOM <- occu(formula = ~ 1
                ~ built + NDVI_median + WVF_PA,
                data = mdata)
coy_best_fit_cov_list <- dredge(coy_SOM, 
                                rank = "AIC")
# plots for deer single occupancy
dem_med_data <- data.frame(DEM_median=seq(60, 180, by=1))
occ.prob.dem_med <- predict(fit_DEM_median, type="state", newdata=dem_med_data, appendData=TRUE)

dem_mdn_plot <- plot(Predicted ~ DEM_median, occ.prob.dem_med, type="l", ylim=c(0,1),
                          xlab="DEM Median",
                          ylab="Expected occupancy probability")
lines(lower ~ DEM_median, occ.prob.dem_med, type="l", col=gray(0.5))
lines(upper ~ DEM_median, occ.prob.dem_med, type="l", col=gray(0.5))

ndvi_mean_data <- data.frame(NDVI_mean=seq(0.15, 0.5, by=0.01))
occ.prob.ndvi_mean <- predict(fit_NDVI_mean, type="state", newdata=ndvi_mean_data, appendData=TRUE)

ndvi_mean_plot <- plot(Predicted ~ NDVI_mean, occ.prob.ndvi_mean, type="l", ylim=c(0,1),
                            xlab="NDVI Mean",
                            ylab="Expected occupancy probability")
lines(lower ~ NDVI_mean, occ.prob.ndvi_mean, type="l", col=gray(0.5))
lines(upper ~ NDVI_mean, occ.prob.ndvi_mean, type="l", col=gray(0.5))

wvf_dist_data <- data.frame(WVF_dist=seq(0, 2000, by=1))
occ.prob.wvf_dist <- predict(fit_WVF, type="state", newdata=wvf_dist_data, appendData=TRUE)

wvf_dist_plot <- plot(Predicted ~ WVF_dist, occ.prob.wvf_dist, type="l", ylim=c(0,1),
                           xlab="WVF Dist",
                           ylab="Expected occupancy probability")
lines(lower ~ WVF_dist, occ.prob.wvf_dist, type="l", col=gray(0.5))
lines(upper ~ WVF_dist, occ.prob.wvf_dist, type="l", col=gray(0.5))

wvf_pa_data <- data.frame(WVF_PA=seq(0, 0.14, by=0.001))
occ.prob.wvf_pa <- predict(fit_WVF_PA, type="state", newdata=wvf_pa_data, appendData=TRUE)

wvf_pa_plot <- plot(Predicted ~ WVF_PA, occ.prob.wvf_pa, type="l", ylim=c(0,1),
                      xlab="WVF PA",
                      ylab="Expected occupancy probability")
lines(lower ~ WVF_PA, occ.prob.wvf_pa, type="l", col=gray(0.5))
lines(upper ~ WVF_PA, occ.prob.wvf_pa, type="l", col=gray(0.5))

###notes for interpretation
## at 1000 buffer, built was the best model , however not under 2 AIC score from null, just 0.62

## for the 500 and 1000 buffer with deer, the best models were WVF Dist, DEM Median, DEM Mean, NDVI Mean, and NDVI Median
## in that order with AIC values all within at most 1.6 difference

## for 2000 buffer with deer, the best models were NDVI Mean, NDVI Median and WVF PA
## in that order with AIC values all within at most 0.58 difference

## for the plots of the deer occupancy, the x-axis of the plots will be needed to changed according to the range
## of the covariates
