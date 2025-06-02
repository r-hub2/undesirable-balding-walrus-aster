
 library(aster)

 data(echinacea)

 # cut and paste from help(predict.aster)
 vars <- c("ld02", "ld03", "ld04", "fl02", "fl03", "fl04",
     "hdct02", "hdct03", "hdct04")
 redata <- reshape(echinacea, varying = list(vars), direction = "long",
     timevar = "varb", times = as.factor(vars), v.names = "resp")
 redata <- data.frame(redata, root = 1)
 pred <- c(0, 1, 2, 1, 2, 3, 4, 5, 6)
 fam <- c(1, 1, 1, 1, 1, 1, 3, 3, 3)
 hdct <- grepl("hdct", as.character(redata$varb))
 redata <- data.frame(redata, hdct = as.integer(hdct))
 level <- gsub("[0-9]", "", as.character(redata$varb))
 redata <- data.frame(redata, level = as.factor(level))
 aout <- aster(resp ~ varb + level : (nsloc + ewloc) + hdct : pop,
     pred, fam, varb, id, root, data = redata)
 # end of cut and paste from help(predict.aster)

 length(aout$dropped)
 # so would have had a problem with restart in version 1.0-3 and before
 aout <- aster(resp ~ varb + level : (nsloc + ewloc) + hdct : pop,
     pred, fam, varb, id, root, data = redata, parm = aout$coefficients,
     optout = TRUE)
 aout$optout$converged
 aout$optout$iterations == 1

