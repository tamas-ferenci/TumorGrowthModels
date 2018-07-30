library( data.table )
library( nlme )
library( lattice )
library( Hmisc )
library( minpack.lm )
library( nlstools )

RawData <- fread( "TeljesAdatok.csv", dec = "," )

RawData$Date <- as.Date( RawData$Date, format = "%Y.%m.%d" )
RawData$MRI <- as.numeric( RawData$MRI )
RawData <- melt( RawData, id.vars = c( "Date", "Code", "Type" ), variable.name = "MeasurementType",
                 value.name = "Volume" )
RawData <- RawData[,.( Date, Type,
                       DateDay = as.numeric( difftime( Date, min( Date ), units = "days" ) ),
                       MeasurementType, Volume ), .( Code ) ]

RawData <- RawData[ !is.na( RawData$Volume ), ]
RawData$Type <- as.factor( RawData$Type )
levels( RawData$MeasurementType ) <- c( "Caliper-1", "Caliper-2", "Caliper-3", "MRI" )
RawData$Code <- as.factor( RawData$Code )
levels( RawData$Code ) <- c( paste0( "Control #", 1:5 ), paste0( "Treated #", 1:9 ) )
RawData$CodeMeasType <- interaction( RawData$Code, RawData$MeasurementType, sep = ", " )

cairo_pdf( "SizeMeasurements.pdf" )
xyplot( Volume ~ DateDay | Code, data = RawData, groups = MeasurementType, xlab = "Time [day]",
        ylab = expression( paste( "Tumor size [ ", mm^3, " ]" ) ), type = "b", layout = c( 7, 2 ),
        auto.key = list( columns = 2, points = FALSE, lines = TRUE ), as.table = TRUE )
dev.off()

mysettings <- list( superpose.line = list( col = rep( trellis.par.get()$superpose.line$col[ 1:4 ],
                                                      each = 2 ),
                                           lty = rep( c( "solid", "dashed"), 4 ) ) )

predgrid <- expand.grid( DateDay = 0:150, Type = unique( RawData$Type ),
                         MeasurementType = unique( RawData$MeasurementType ) )


# Exponential

## I use this approach instead of nlsList so that I can use nlsLM
IndExp <- do.call( rbind, lapply( unique( RawData$CodeMeasType ), function( cmt ) {
  fit <- nlsLM( Volume ~ V0*exp( a*DateDay ), data = RawData[ CodeMeasType==cmt, ],
                start = c( V0 = 20, a = 0.2 ),
                control = nls.lm.control( maxiter = 1024, maxfev = 10000 ) )
  data.table( cmt = cmt, par = names( coef( fit ) ), est = coef( fit ), confint2( fit ) )
} ) )

cairo_pdf( "IndExp.pdf" )
Dotplot( as.numeric( cmt ) ~ Cbind( est, `2.5%`, `97.5 %` ) | par, data = IndExp,
         scales = list( relation = "free",
                        y = list( at = 1:length( unique( IndExp$cmt ) ),
                                  labels = levels( IndExp$cmt )[ 1:length( unique( IndExp$cmt ) ) ]
                        ) ), xlab = "", ylab = "" )
dev.off()

fit <- nlme( model = Volume ~ V0*exp( a*DateDay ), data = RawData,
             fixed = V0 + a ~ Type + MeasurementType, random = list( Code = pdDiag( V0 + a ~ 1 ) ),
             start = c( rep( 1, 5 ), rep( 0.1, 5 ) ), 
             ## these require more or less 'trial and error'...
             control = nlmeControl( maxIter = 100, msMaxIter = 200 ) )

cairo_pdf( "ExpHeterosced.pdf" )
plot( fit, resid( ., type = "normalized" ) ~ log( fitted(.) )|MeasurementType )
dev.off()

fit <- update( fit, weights = varPower( form = ~fitted(.)|MeasurementType ) )

cairo_pdf( "ExpHeterosced2.pdf" )
plot( fit, resid( ., type = "normalized" ) ~ log( fitted(.) )|MeasurementType )
dev.off()

cairo_pdf( "ExpAutocorr.pdf" )
plot( ACF( fit, resType = "normalized", maxLag = 50 ), alpha = 0.05 )
dev.off()

fit <- update( fit, corr = corARMA( form = ~DateDay|Code/MeasurementType, p = 1, q = 2 ) )

cairo_pdf( "ExpAutocorr2.pdf" )
plot( ACF( fit, resType = "normalized", maxLag = 50 ), alpha = 0.05 )
dev.off()

summary( fit )

cairo_pdf( "PredExponential.pdf" )
xyplot( pred ~ DateDay, data = cbind( predgrid, pred = predict( fit, predgrid, level = 0 ) ),
        xlab = "Time [day]", ylab = expression( paste( "Tumor size [ ", mm^3, " ]" ) ),
        groups = interaction( Type, MeasurementType ), type = "l",
        auto.key = list( text = as.character( rep( unique( RawData$MeasurementType ), each = 2 ) ),
                         columns = 4, lines = TRUE, points = FALSE ),
        par.settings = mysettings )
dev.off()

# Gompertz

IndGompertz <- do.call( rbind, lapply( unique( RawData$CodeMeasType ), function( cmt ) {
  fit <- nlsLM( Volume ~ SSgompertz( DateDay, Asym, b2, b3 ), data = RawData[ CodeMeasType==cmt, ],
                start = c( Asym = 50000, b2 = 10, b3 = 0.8 ),
                control = nls.lm.control( maxiter = 1024, maxfev = 10000 ) )
  data.table( cmt = cmt, par = names( coef( fit ) ), est = coef( fit ),
              se = sqrt( diag( vcov( fit ) ) ), confint2( fit ) )
} ) )

cairo_pdf( "IndGompertz.pdf" )
Dotplot( as.numeric( cmt ) ~ Cbind( est, `2.5 %`, `97.5 %` ) | par, data = IndGompertz,
         scales = list( relation = "free",
                        y = list( at = 1:length( unique( IndGompertz$cmt ) ),
                                  labels = levels( IndGompertz$cmt )[ 1:length(
                                    unique( IndGompertz$cmt ) ) ]
                        ) ), xlab = "", ylab = "", layout = c( 3, 1 ) )
dev.off()

histogram( ~se | par, data = IndGompertz, scales = list( relation = "free",x = list( log = 10 ) ) )
IndGompertz[ par=="Asym"&se>1e12 ]
RawDataGompertz <- RawData[ !CodeMeasType%in%IndGompertz[ (par=="Asym"&se>1e10)|(par=="b2"&se>300) ]$cmt ]
RawDataGompertz$CodeMeasType <- droplevels( RawDataGompertz$CodeMeasType )
## While it is true that the mixed-effects models make use of all data, even those that don't fit
## to the same curve, these are so extremely outlying observations that we're probably better off
## without them (they'd make the model estimation, which is already a tough task from so
## early-phase data, even more problematic).

IndGompertz <- do.call( rbind, lapply( unique( RawDataGompertz$CodeMeasType ), function( cmt ) {
  fit <- nlsLM( Volume ~ SSgompertz( DateDay, Asym, b2, b3 ),
                data = RawDataGompertz[ CodeMeasType==cmt, ],
                start = c( Asym = 50000, b2 = 10, b3 = 0.8 ),
                control = nls.lm.control( maxiter = 1024, maxfev = 10000 ) )
  data.table( cmt = cmt, par = names( coef( fit ) ), est = coef( fit ),
              se = sqrt( diag( vcov( fit ) ) ), confint2( fit ) )
} ) )

cairo_pdf( "IndGompertz2.pdf" )
Dotplot( as.numeric( cmt ) ~ Cbind( est, `2.5 %`, `97.5 %` ) | par, data = IndGompertz,
         scales = list( relation = "free",
                        y = list( at = 1:length( unique( IndGompertz$cmt ) ),
                                  labels = levels( IndGompertz$cmt )[ 1:length(
                                    unique( IndGompertz$cmt ) ) ]
                        ) ), xlab = "", ylab = "", layout = c( 3, 1 ) )
dev.off()

fit <- nlme( model = Volume ~ SSgompertz( DateDay, Asym, b2, b3 ), data = RawDataGompertz,
             fixed = Asym + b2 + b3 ~ Type + MeasurementType,
             random = list( Code = pdDiag( Asym + b2 + b3 ~ 1 ) ),
             start = c( 100000, rep( 0.1, 4 ), 10, rep( 0.1, 4 ), 0, rep( 0.1, 4 ) ),
             control = nlmeControl( maxIter = 100, msMaxIter = 100 ) )

cairo_pdf( "GompertzHeterosced.pdf" )
plot( fit, resid( ., type = "normalized" ) ~ log( fitted(.) )|MeasurementType )
dev.off()

fit <- update( fit, weights = varPower( form = ~fitted(.)|MeasurementType ) )

cairo_pdf( "GompertzHeterosced2.pdf" )
plot( fit, resid( ., type = "normalized" ) ~ log( fitted(.) )|MeasurementType )
dev.off()

cairo_pdf( "GompertzAutocorr.pdf" )
plot( ACF( fit, resType = "normalized", maxLag = 50 ), alpha = 0.05 )
dev.off()

fit <- update( fit, corr = corARMA( form = ~DateDay|Code/MeasurementType, p = 1, q = 2 ) )

cairo_pdf( "GompertzAutocorr2.pdf" )
plot( ACF( fit, resType = "normalized", maxLag = 50 ), alpha = 0.05 )
dev.off()

summary( fit )

cairo_pdf( "PredGompertz.pdf" )
xyplot( pred ~ DateDay, data = cbind( predgrid, pred = predict( fit, predgrid, level = 0 ) ),
        xlab = "Time [day]", ylab = expression( paste( "Tumor size [ ", mm^3, " ]" ) ),
        groups = interaction( Type, MeasurementType ), type = "l",
        auto.key = list( text = as.character( rep( unique( RawData$MeasurementType ), each = 2 ) ),
                         columns = 4, lines = TRUE, points = FALSE ),
        par.settings = mysettings )
dev.off()

# Logistic

IndLogistic <- do.call( rbind, lapply( unique( RawData$CodeMeasType ), function( cmt ) {
  fit <- nlsLM( Volume ~ SSlogis( DateDay, Asym, xmid, scal ), data = RawData[ CodeMeasType==cmt, ],
                start = c( Asym = 20000, xmid = 15, scal = 3 ),
                control = nls.lm.control( maxiter = 1024, maxfev = 10000 ) )
  data.table( cmt = cmt, par = names( coef( fit ) ), est = coef( fit ),
              se = sqrt( diag( vcov( fit ) ) ), confint2( fit ) )
} ) )

cairo_pdf( "IndLogistic.pdf" )
Dotplot( as.numeric( cmt ) ~ Cbind( est, `2.5 %`, `97.5 %` ) | par, data = IndLogistic,
         scales = list( relation = "free",
                        y = list( at = 1:length( unique( IndLogistic$cmt ) ),
                                  labels = levels( IndLogistic$cmt )[ 1:length(
                                    unique( IndLogistic$cmt ) ) ]
                        ) ), xlab = "", ylab = "", layout = c( 3, 1 ) )
dev.off()

histogram( ~se | par, data = IndLogistic, scales = list( relation = "free",x = list( log = 10 ) ) )
IndLogistic[ par=="Asym"&se>1e7 ]
RawDataLogistic <- RawData[ !CodeMeasType%in%IndLogistic[ par=="Asym"&se>1e7 ]$cmt ]
RawDataLogistic$CodeMeasType <- droplevels( RawDataLogistic$CodeMeasType )

IndLogistic <- do.call( rbind, lapply( unique( RawDataLogistic$CodeMeasType ), function( cmt ) {
  fit <- nlsLM( Volume ~ SSlogis( DateDay, Asym, xmid, scal ),
                data = RawDataLogistic[ CodeMeasType==cmt, ],
                start = c( Asym = 20000, xmid = 15, scal = 3 ),
                control = nls.lm.control( maxiter = 1024, maxfev = 10000 ) )
  data.table( cmt = cmt, par = names( coef( fit ) ), est = coef( fit ),
              se = sqrt( diag( vcov( fit ) ) ), confint2( fit ) )
} ) )

cairo_pdf( "IndLogistic2.pdf" )
Dotplot( as.numeric( cmt ) ~ Cbind( est, `2.5 %`, `97.5 %` ) | par, data = IndLogistic,
         scales = list( relation = "free",
                        y = list( at = 1:length( unique( IndLogistic$cmt ) ),
                                  labels = levels( IndLogistic$cmt )[ 1:length(
                                    unique( IndLogistic$cmt ) ) ]
                        ) ), xlab = "", ylab = "", layout = c( 3, 1 ) )
dev.off()

fit <- nlme( model = Volume ~ SSlogis( DateDay, Asym, xmid, scal ), data = RawDataLogistic,
             fixed = Asym + xmid + scal ~ Type + MeasurementType,
             random = list( Code = pdDiag( Asym + xmid + scal ~ 1 ) ),
             start = c( 20000, rep( 0, 4 ), 15, rep( 0, 4 ), 3, rep( 0, 4 ) ),
             control = nlmeControl( maxIter = 100, msMaxIter = 100 ) )

cairo_pdf( "LogisticHeterosced.pdf" )
plot( fit, resid( ., type = "normalized" ) ~ log( fitted(.) )|MeasurementType )
dev.off()

fit <- update( fit, weights = varPower( form = ~fitted(.)|MeasurementType ) )

cairo_pdf( "LogisticHeterosced2.pdf" )
plot( fit, resid( ., type = "normalized" ) ~ log( fitted(.) )|MeasurementType )
dev.off()

cairo_pdf( "LogisticAutocorr.pdf" )
plot( ACF( fit, resType = "normalized", maxLag = 50 ), alpha = 0.05 )
dev.off()

fit <- update( fit, corr = corARMA( form = ~DateDay|Code/MeasurementType, p = 1, q = 2 ) )

cairo_pdf( "LogisticAutocorr2.pdf" )
plot( ACF( fit, resType = "normalized", maxLag = 50 ), alpha = 0.05 )
dev.off()

summary( fit )

cairo_pdf( "PredLogistic.pdf" )
xyplot( pred ~ DateDay, data = cbind( predgrid, pred = predict( fit, predgrid, level = 0 ) ),
        xlab = "Time [day]", ylab = expression( paste( "Tumor size [ ", mm^3, " ]" ) ),
        groups = interaction( Type, MeasurementType ), type = "l",
        auto.key = list( text = as.character( rep( unique( RawData$MeasurementType ), each = 2 ) ),
                         columns = 4, lines = TRUE, points = FALSE ),
        par.settings = mysettings )
dev.off()
