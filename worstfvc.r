library(sclero)
library(plyr)
library(ggplot2)

data(tPFT)
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

tPFT <- transform(tPFT, days=dtdiff(first_seen(PtID), Date))
tPFT <- transform(tPFT, years=days %/% 365)

worstobs <- function(xs) {
  xs <- na.omit(xs)
  ok <- length(xs) > 0
  result <- if (ok) median(quantile(xs, 0.1)) else NA
  return(result)
}

yearlb <- 4
targets <- ddply(subset(tPFT, years >= yearlb), ~ PtID, summarize,
                 lowpfvc=worstobs(perc.FVC.of.predicted))

patients <- targets$PtID
pft <- subset(tPFT, PtID %in% patients & years < yearlb)

yearavgs <- function(df, varname, maxyr) {
  s <- rep(0, maxyr)
  c <- rep(0, maxyr)
  
  for (i in seq_along(nrow(df))) {
    y <- df$years[i]
    v <- df[[varname]][i]
    if (is.na(v)) next
    s[y+1] <- s[y+1] + df[[varname]][i]
    c[y+1] <- c[y+1] + 1
  }
  
  df <- data.frame(year=(1:maxyr - 1))
  s[c > 0] <- s[c > 0] / c[c > 0]
  df[[varname]] <- s
  
  return(df)
}

yearly <- ddply(pft, ~ PtID, yearavgs, "perc.FVC.of.predicted", yearlb)

require(reshape2)

yearly <- dcast(yearly, PtID ~ year, value.var="perc.FVC.of.predicted")
names(yearly) <- c("PtID", paste0("Y", names(yearly)[-1]))

addindicator <- function(df) {
  nms <- paste0("obs.", names(yearly)[-1])
  obs <- data.frame(df[, 2] > 0, df[, 3] > 0, df[, 4] > 0, df[, 5] > 0)
  names(obs) <- nms
  df <- cbind(df, obs)
  return(df)
}

yearly <- ddply(yearly, ~ PtID, addindicator)
yearly <- merge(yearly, targets, by="PtID", all.x=TRUE)

yearly <- yearly[yearly$obs.Y0, ]

m0 <- lm(lowpfvc ~ Y0, data=yearly)
m1 <- lm(lowpfvc ~ Y0 + Y1, data=yearly)
m2 <- lm(lowpfvc ~ Y0 + Y1 + Y2, data=yearly)
m3 <- lm(lowpfvc ~ Y0 + Y1 + Y2 + Y3, data=yearly)

yearly <- transform(yearly, n.obs = obs.Y1 + obs.Y2 + obs.Y3)
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

residuals <- yearly[, c("PtID", "lowpfvc", "m0", "m1", "m2", "m3")]
residuals <- melt(residuals, id.vars=c("PtID", "lowpfvc"),
                  variable.name="model", value.name="prediction")
res2 <- melt(residuals, id.vars=c("PtID", "lowpfvc", "m0"),
                  variable.name="model", value.name="prediction")
residuals[["group"]] <- factor(grplbl(residuals$PtID), 0:3)
res2[["group"]] <- factor(grplbl(residuals$PtID), 0:3)
res2 <- transform(res2, baseline=lowpfvc - m0)
res2 <- transform(res2, extended=lowpfvc - prediction)

p <- ggplot(res2, aes(x = abs(baseline), y = abs(extended)))
p <- p + geom_point(aes(color=group)) + facet_wrap(~ model, ncol=2)
p + labs(title="Residuals of Baseline against Extended Models")

p <- ggplot(residuals, aes(x = group == 0, y = lowpfvc - prediction, color=group))
p <- p + geom_jitter() + facet_wrap(~ model, ncol=2)
p + labs(title="Residuals Faceted on Model (0-3)")

group3 <- subset(residuals, group == 3)
group3 <- arrange(group3, PtID)
group3[["years"]] <- 0:3
group3 <- transform(group3, PtID = as.factor(PtID))
above15 <- subset(group3, model=="m0" & abs(lowpfvc-prediction) > 15)$PtID
p <- ggplot(subset(group3, PtID %in% above15), aes(x=years, y=abs(lowpfvc-prediction)))
p <- p + geom_line(aes(color=PtID), alpha=0.3) + geom_smooth(method=lm)
p + labs(title="Years of Data against Absolute Error")
