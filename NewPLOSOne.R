library( data.table )
library( lattice )
library( nlme )

RawData <- rbind(
  melt( cbind( Time = as.numeric( fread( "timeERC_PI.csv" )$V1 ), fread( "dataERC_PI.csv" ),
               Study = "HT-29", Treatment = "metronomic" ), id.vars = c( "Study", "Treatment", "Time" ) ),
  melt( cbind( Time = as.numeric( fread( "timeERC_PIC.csv" )$V1 ), fread( "dataERC_PIC.csv" ),
               Study = "HT-29", Treatment = "no treatment" ), id.vars = c( "Study", "Treatment", "Time" ) ),
  melt( cbind( Time = as.numeric( fread( "timePI.csv" )$V1 ), fread( "dataPI.csv" ),
               Study = "C38", Treatment = "no treatment" ), id.vars = c( "Study", "Treatment", "Time" ) ),
  melt( cbind( Time = as.numeric( fread( "timePIII_2.csv" )$V1 ), fread( "dataPIII_2.csv" ),
               Study = "C38", Treatment = "metronomic" ), id.vars = c( "Study", "Treatment", "Time" ) ) )

names( RawData )[ names( RawData )=="variable" ] <- "Mouse"
RawData$Mouse <- paste0( RawData$Study, ", ", RawData$Treatment, ", ", substring( RawData$Mouse, 2 ) )
RawData <- RawData[ !is.nan( value ) ]
RawData$logvolume <- log( RawData$value )
RawData <- RawData[ order( Mouse, Time ) ]
RawData$First <- ifelse( duplicated( RawData$Mouse ), "Notfirst", "First" )
RawData <- RawData[ order( Mouse, -Time ) ]
RawData$Last <- ifelse( duplicated( RawData$Mouse ), "Notlast", "Last" )
RawData <- RawData[ order( Mouse, Time ) ]
RawData$Treatment <- relevel( as.factor( RawData$Treatment ), ref = "no treatment" )

cairo_pdf( "Baseline.pdf" )
dotplot( logvolume ~ Treatment | Study, data = RawData[First=="First"], ylab = "Log volume" )
dev.off()

coin::wilcox_test( logvolume ~ Treatment, data = RawData[Study=="C38"&First=="First"], distribution = "exact" )
coin::wilcox_test( logvolume ~ Treatment, data = RawData[Study=="HT-29"&First=="First"], distribution = "exact" )

cairo_pdf( "GrowthCurves.pdf" )
xyplot( logvolume ~ Time | Treatment+Study, groups = Mouse, data = RawData, type = "l",
        ylab = "Log volume", xlab = "Time [days]" )
dev.off()

fit_C38_1 <- summary( lm( Notfirst ~ Treatment, data = dcast( RawData[Study=="C38"&(First=="First"|Last=="Last")],
                                                              Treatment+Study+Mouse ~ First, value.var = "logvolume" ) ) )
fit_HT29_1 <- summary( lm( Notfirst ~ Treatment, data = dcast( RawData[Study=="HT-29"&(First=="First"|Last=="Last")],
                                                               Treatment+Study+Mouse ~ First, value.var = "logvolume" ) ) )

fit_C38_2 <- summary( lm( Notfirst ~ Treatment+First, data = dcast( RawData[Study=="C38"&(First=="First"|Last=="Last")],
                                                                    Treatment+Study+Mouse ~ First, value.var = "logvolume" ) ) )
fit_HT29_2 <- summary( lm( Notfirst ~ Treatment+First, data = dcast( RawData[Study=="HT-29"&(First=="First"|Last=="Last")],
                                                                     Treatment+Study+Mouse ~ First, value.var = "logvolume" ) ) )

fit_C38_3 <- lme( logvolume ~ Time*Treatment, data = RawData[Study=="C38"], random = ~Time|Mouse )
fit_HT29_3 <- lme( logvolume ~ Time*Treatment, data = RawData[Study=="HT-29"], random = ~Time|Mouse )


res <- cbind(
  rbind(
    fit_C38_1$coefficients[ "Treatmentmetronomic", c( "Estimate", "Std. Error", "Pr(>|t|)" ) ],
    fit_C38_2$coefficients[ "Treatmentmetronomic", c( "Estimate", "Std. Error", "Pr(>|t|)" ) ],
    summary(fit_C38_3)$tTable[ "Time:Treatmentmetronomic", c( "Value", "Std.Error", "p-value" ) ]
  ),
  rbind(
    fit_HT29_1$coefficients[ "Treatmentmetronomic", c( "Estimate", "Std. Error", "Pr(>|t|)" ) ],
    fit_HT29_2$coefficients[ "Treatmentmetronomic", c( "Estimate", "Std. Error", "Pr(>|t|)" ) ],
    summary(fit_HT29_3)$tTable[ "Time:Treatmentmetronomic", c( "Value", "Std.Error", "p-value" ) ]
  )
)

apply( round(res,4), 1, function(x) paste(x,collapse="&") )

res2 <- data.frame(
  cbind( intervals(fit_C38_3)$fixed[,c("est.","lower","upper")], p=summary(fit_C38_3)$tTable[,"p-value"] ),
  cbind( intervals(fit_HT29_3)$fixed[,c("est.","lower","upper")], p=summary(fit_HT29_3)$tTable[,"p-value"] )
)
res2[,4] <- Hmisc::format.pval(res2[,4],eps=0.001)
res2[,8] <- Hmisc::format.pval(res2[,8],eps=0.001)
res2[,c(1:3,5:7)] <- round(res2[,c(1:3,5:7)],4)

apply( res2, 1, function(x) paste(x,collapse="&") )

summary(fit_C38_3)
summary(fit_HT29_3)

log(2)/(intervals(fit_C38_3)$fixed["Time","est."]+
          c(0,intervals(fit_C38_3)$fixed["Time:Treatmentmetronomic",c("est.","lower","upper")]))
log(2)/(intervals(fit_HT29_3)$fixed["Time","est."]+
          c(0,intervals(fit_HT29_3)$fixed["Time:Treatmentmetronomic",c("est.","lower","upper")]))

fit_full <- lme( logvolume ~ Time*Treatment*Study, data = RawData, random = ~Time|Mouse )
summary(fit_full)
fit_full_2 <- update( fit_full, fixed = .~.-Time:Treatment:Study )
summary(fit_full_2)
intervals(fit_full_2)