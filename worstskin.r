library(sclero)
library(plyr)
library(ggplot2)
library(reshape2)

data(tVisit)
data(tPtData)

first_seen <- local({
  pids <- tPtData$PtID
  dates <- tPtData$DateFirstSeen
  lookup <- structure(dates, names=pids)
  function(pid) lookup[as.character(pid)]
})

dtdiff <- function(from, to) {
  as.integer(as.Date(to) - as.Date(from))
}

tVisit <- transform(tVisit, days=dtdiff(first_seen(PtID), Visit.Date))
tVisit <- transform(tVisit, years=days %/% 365)

worstobs <- function(xs) {
  xs <- na.omit(xs)
  ok <- length(xs) > 0
  result <- if (ok) median(quantile(xs, 0.9)) else NA
  return(result)
}

yearlb <- 4
targets <- ddply(subset(tVisit, years >= yearlb), ~ PtID, summarize,
                 highskin=worstobs(Total.Skin.Score))

patients <- targets$PtID
clinic <- subset(tVisit, PtID %in% patients & years < yearlb)
clinic <- clinic[, c("PtID", "years", "Total.Skin.Score")]
clinic$years[clinic$years < 0] <- 0

yearly <- ddply(na.omit(clinic), ~ PtID + years, summarize,
                Total.Skin.Score = mean(Total.Skin.Score))
yearly <- dcast(yearly, PtID ~ years, value.var="Total.Skin.Score")
names(yearly) <- c("PtID", paste0("Y", names(yearly)[-1]))

addindicator <- function(df) {
  nms <- paste0("obs.", names(yearly)[-1])
  obs <- data.frame(!is.na(df[, 2]), !is.na(df[, 3]), !is.na(df[, 4]), !is.na(df[, 5]))
  names(obs) <- nms
  df <- cbind(df, obs)
  return(df)
}

yearly <- ddply(yearly, ~ PtID, addindicator)
yearly <- merge(yearly, targets, by="PtID", all.x=TRUE)

yearly <- yearly[yearly$obs.Y0, ]
yearly$Y1[is.na(yearly$Y1)] <- 0
yearly$Y2[is.na(yearly$Y2)] <- 0
yearly$Y3[is.na(yearly$Y3)] <- 0

m0 <- lm(highskin ~ Y0, data=yearly)
m1 <- lm(highskin ~ Y0 + I(Y1*obs.Y1), data=yearly)
m2 <- lm(highskin ~ Y0 + I(Y1*obs.Y1) + I(Y2*obs.Y2), data=yearly)
m3 <- lm(highskin ~ Y0 + I(Y1*obs.Y1) + I(Y2*obs.Y2) + I(Y3*obs.Y3), data=yearly)

p1 <- subset(yearly, obs.Y1)$PtID
p2 <- subset(yearly, obs.Y1 & obs.Y2)$PtID
p3 <- subset(yearly, obs.Y1 & obs.Y2 & obs.Y3)$PtID

loocv <- function(i, df, m) {
  fit <- lm(formula(m), data=df[-i, ])
  predict(fit, df[i, ])
}

idx <- 1:nrow(yearly)
yearly[["m0"]] <- sapply(idx, function(i) loocv(i, yearly, m0))
yearly[["m1"]] <- sapply(idx, function(i) loocv(i, yearly, m1))
yearly[["m2"]] <- sapply(idx, function(i) loocv(i, yearly, m2))
yearly[["m3"]] <- sapply(idx, function(i) loocv(i, yearly, m3))

grplbl <- function(pid) {
  ifelse(pid %in% p3, 3,
         ifelse(pid %in% p2, 2,
                ifelse(pid %in% p1, 1, 0)))
}

residuals <- yearly[, c("PtID", "highskin", "m0", "m1", "m2", "m3")]
residuals <- melt(residuals, id.vars=c("PtID", "highskin"),
                  variable.name="model", value.name="prediction")
res2 <- melt(residuals, id.vars=c("PtID", "highskin", "m0"),
                  variable.name="model", value.name="prediction")
residuals[["group"]] <- factor(grplbl(residuals$PtID), 0:3)
res2[["group"]] <- factor(grplbl(residuals$PtID), 0:3)
res2 <- transform(res2, baseline=highskin - m0)
res2 <- transform(res2, extended=highskin - prediction)

p <- ggplot(res2, aes(x = abs(baseline), y = abs(extended)))
p <- p + geom_point(aes(color=group)) + facet_wrap(~ model, ncol=2)
p + labs(title="Residuals of Baseline against Extended Models")

p <- ggplot(residuals, aes(x = group == 0, y = highskin - prediction, color=group))
p <- p + geom_jitter() + facet_wrap(~ model, ncol=2)
p + labs(title="Residuals Faceted on Model (0-3)")

group3 <- subset(residuals, group == 3)
group3 <- arrange(group3, PtID)
group3[["years"]] <- 0:3
group3 <- transform(group3, PtID = as.factor(PtID))
above15 <- subset(group3, model=="m0" & abs(highskin-prediction) > 15)$PtID
p <- ggplot(subset(group3, PtID %in% above15), aes(x=years, y=abs(highskin-prediction)))
p <- p + geom_line(aes(color=PtID), alpha=0.3) + geom_smooth(method=lm)
p + labs(title="Years of Data against Absolute Error")
