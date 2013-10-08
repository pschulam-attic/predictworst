library(sclero)
library(ggplot2)

data(tPtData)
data(sclero.combined)
data(lifetime.worst)

sclero.combined <- merge(sclero.combined, tPtData, by="PtID", all.x=TRUE)

get_baseline <- function(df, max.days) {
  require(plyr)

  df <- arrange(df, Date)
  first.seen <- as.Date(df$DateFirstSeen)[1]
  
  days.since <- as.integer(as.Date(df$Date) - first.seen)
  base.df <- df[days.since <= max.days, ]
  
  var.names <- c("Height", "Weight", "age", "Sex", "Race1",
                 "perc.FVC.of.predicted", "perc.DLCO.of.predicted",
                 "Total.Skin.Score", "GI.Sev.Score", "RP.Sev.Score",
                 "ANA", "Scl.70", "ACA")
  
  base.df[, var.names]
}

library(plyr)

base.visits <- ddply(sclero.combined, ~ PtID, get_baseline, 1.5 * 365)
base.visits <- transform(base.visits,
                         PtID = as.factor(PtID),
                         ANA = factor(ANA, 0:1),
                         Scl.70 = factor(Scl.70, 0:1),
                         ACA = factor(ACA, 0:1),
                         Sex = factor(Sex, 0:1),
                         Afr = factor(ifelse(Race1 == 2, 1, 0), 0:1))

maxna <- function(xs) {
  xs <- na.omit(xs)
  if (length(xs) > 0) {
    if (is.factor(xs)) {
      result <- xs[length(xs)]
    } else {
      result <- mean(xs)
    }
  } else {
    result <- NA
  }
  return(result)
}

summarize_baseline <- function(df) {
  require(plyr)
  
  result <- lapply(names(df), function(n) maxna(df[[n]]))
  names(result) <- names(df)
  result <- as.data.frame(result)
  return(result)
}

pred.data <- ddply(base.visits, ~ PtID, summarize_baseline)
pred.data <- merge(pred.data, lifetime.worst, by="PtID", all.x=TRUE)

loocv <- function(i, df, glmfit, fam) {
  fit <- glm(formula(glmfit), data=df[-i, ], family=fam)
  p <- predict(fit, df[i, ], se.fit=TRUE, type="response")
  c(p$fit, p$se.fit, p$residual.scale, df$will.decline[i])
}

df <- na.omit(pred.data)

df <- transform(
  df, will.decline=ifelse(perc.FVC.of.predicted > 70 & worst <= 70, 1, 0))

healthy <- subset(df, perc.FVC.of.predicted > 70)
binfit <- glm(will.decline ~ Height + Weight + age + Sex +
                Afr + perc.FVC.of.predicted + 
                perc.DLCO.of.predicted + Total.Skin.Score +
                GI.Sev.Score + RP.Sev.Score + ACA + Scl.70, family="binomial",
              data=healthy)

binfit.cv <- sapply(1:nrow(healthy), function(i) loocv(i, healthy, binfit, "binomial"))

thresh_predict <- function(t, p) {
  dec <- ifelse(p >= t, 1, 0)
  return(dec)
}

det_curve <- function(y, y.hat) {
  fa <- sum(y.hat[y == 0]) / sum(y == 0)
  md <- sum(1 - y.hat[y == 1]) / sum(y == 1)
  c(fa=fa, md=md)
}

roc_curve <- function(y, y.hat) {
  tp <- sum(y)
  tn <- sum(1 - y)
  tpr <- sum(y.hat[y == 1]) / tp
  fpr <- sum(y.hat[y == 0]) / tn
  c(fpr=fpr, tpr=tpr)
}

score_binfit <- function(t) {
  y <- binfit.cv[4, ]
  p <- binfit.cv[1, ]
  d <- thresh_predict(t, p)
  s <- roc_curve(y, d)
  return(s)
}

ts <- seq(0, 1, length.out=1000)
scores <- sapply(ts, score_binfit)
p <- qplot(scores[1, ], scores[2, ], xlab="False Positive Rate",
           ylab="True Positive Rate")
p + labs(title="ROC Curve for Identifying At-risk Patients")

ranks <- order(abs(scores[1, ] - scores[2, ]))
