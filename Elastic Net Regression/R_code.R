#import required libraries
rm(list=ls())
library(caret)
library(glmnet)
library(MASS)
library(mlbench)
library(psych)


# loading the database and checking the structure
data <- Boston
# The structure and some graphical summaries of the data
head(data)
dim(data)
str(data)
pairs.panels(data[c(-4,-14)],cex=2)  #column 4 is a dummy variable and 14 is the response.

# Train-Test split
# set seed
set.seed(222)
n <- nrow(data)
# Generate 1 and 2 randomly for 506 times, and the probability of 1 is 0.7, 
# the probability of 1 is 0.3.
ind <- sample(2, n, replace = TRUE, prob = c(0.7, 0.3))
train <- data[ind==1,] # training date
test <- data[ind==2,] # test data
dim(train)
dim(test)

# Custom Control Parameters
custom <- trainControl(method = "cv", number = 10, verboseIter = TRUE)

# Linear Model
set.seed(1234)
lm<- train(medv~., train, method='lm', trControl=custom)
lm$results
lm
summary(lm)
par(mfrow=c(2,2),mar=c(2,2,2,2))
plot(lm$finalModel)

# Ridge Rrgression
set.seed(1234)
ridge <- train(medv~.,
               train,
               method = 'glmnet',
               tuneGrid = expand.grid(alpha = 0,
                                      lambda = seq(0, 1, length = 100)),
               trControl = custom)
ridge$results
plot(ridge)
windows(width = 10, height = 6)
plot(ridge$finalModel, xvar = "lambda", label = T)
plot(ridge$finalModel, xvar = 'dev', label=T)
plot(varImp(ridge, scale=T))

# Lasso Rrgression
set.seed(1234)
lasso <- train(medv~.,
               train,
               method = 'glmnet',
               tuneGrid = expand.grid(alpha = 1,
                                      lambda = seq(0.0001, 0.1, length = 100)),
               trControl = custom)

lasso$results
plot(lasso)
windows(width = 10, height = 6)
plot(lasso$finalModel, xvar = 'lambda', label=T)
plot(lasso$finalModel, xvar = 'dev', label=T)
plot(varImp(lasso,scale = T))

# Elastic Net
set.seed(1234)
en <- train(medv~.,
            train,
            method = 'glmnet',
            tuneGrid = expand.grid(alpha = seq(0,1,length=10),
                                   lambda = seq(0,0.2,length=10)),
            trControl = custom)

en$results
en_table <- data.frame(en$results)
en_table$lambda <- factor(round(en_table[,2],3))
ggplot(en_table, aes(x=alpha, y=RMSE,color=lambda)) +
  geom_point() +
  geom_line(size=0.5) +                        
  scale_colour_manual(values=c("red","blue","darkgreen","gray","orange","purple","yellow","green", "darkblue","darkred"))+
  labs(title="Regularization Parameter", x = "Mixing Percentage", y = "RMSE (Cross-Validation)") +
  theme_classic()+
  theme(legend.position = 'top',plot.title = element_text(hjust=0.5))
windows(width = 10, height = 6)
plot(en$finalModel, xvar = 'lambda', label=T)
plot(en$finalModel, xvar = 'dev', label=T)
plot(varImp(en))


# Save Models
saveRDS(lm, "linear_regression.rds")
mlin <- readRDS("linear_regression.rds")
saveRDS(ridge, "ridge.rds")
mrid <- readRDS("ridge.rds")
saveRDS(lasso, "lasso.rds")
mlas <- readRDS("lasso.rds")
saveRDS(en, "en.rds")
men <- readRDS("en.rds")

# Compare models using MSE
RMSE_train <- data.frame(matrix(nrow=1000, ncol=4))
mlist <- list(mlin, mrid, mlas, men)
for (i in 1:1000){
  set.seed(i)
  boost_train <- sample(train,nrow(train), replace = TRUE)
  for (j in 1:4) {
    p <- predict(mlist[j], boost_train)
    RMSE_train[i,j]<-sqrt(mean((boost_train$medv-p[[1]])^2))
  }
}

Mean_RMSE_train <- c(mean(RMSE_train[,1]),mean(RMSE_train[,2]),mean(RMSE_train[,3]),mean(RMSE_train[,4]))


RMSE_test <- data.frame(matrix(nrow=1000, ncol=4))
mlist <- list(mlin, mrid, mlas, men)
for (i in 1:1000){
  set.seed(i)
  boost_test <- sample(test,nrow(test), replace = TRUE)
  for (j in 1:4) {
    p <- predict(mlist[j], boost_test)
    RMSE_test[i,j]<-sqrt(mean((boost_test$medv-p[[1]])^2))
    }
}

Mean_RMSE_test <- c(mean(RMSE_test[,1]),mean(RMSE_test[,2]),mean(RMSE_test[,3]),mean(RMSE_test[,4]))







