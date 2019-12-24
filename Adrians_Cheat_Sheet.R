#########################################################
##########################
####ADRIAN'S ML CHEAT SHEET
#######Formulas the computer memorizes so you don't have to
#########################
########################################################

############Calculate mean square error for a continuous variable prediction
#Takes in two vectors - y(test), and y(predicted)
mean_square_error <- function(y_test, y_predicted){
  
  mse <- ((sum(y_predicted - y_test))^2)/nrow(y_test)  
  
  return(paste("MSE: ", mse, " RMSE: ", mse^0.5))
  
}

##

############Create a confusion matrix
#Takes in two vectors - y(test), and y(predicted)
confusion_matrix <- function(y_test, y_predicted){
  
  confusion_matrix <- table(y_test, y-predicted)
  
  return(confusion_matrix)
  
}

############Create an elbow plot for a KNN
#Specify your max_k, and split your training/testing data into x and y. training_y should be a vector.
elbow_plot <- function(max_k, training_x, testing_x, training_y){
  predicted_y <- NULL
  error_rate <- NULL

  for (i in 1:max_k){
    predicted_y[i] <- knn(training_x, testing_x, training_y, k=i) 
    error_rate[i] <- mean(predicted_y != training_y)
  }
  
  elbow_plot <- data.frame(cbind(k=c(1:max_k),error_rate))
  return(ggplot(data=elbow_plot, aes(x=k, y=error_rate)) +
         geom_point() +
         geom_line(lty="dotted"))
  
}

############Calculate misclassification rate and accuracy for a logistic regression 
accuracy <- function(y_predicted, y_test){
  misClasificError <- mean(y_predicted != y_test)
  return(paste("Misclassification Error: ", misClasificError, ' Accuracy: ',1-misClasificError))
  
}

###########Auto_scaler
####Uses the maxima and minima of columns to scale them accordingly
#Argument is a dataframe
auto_scaler <- function(dataframe){
  column_maxes <- apply(dataframe, 2, max)
  column_mins <- apply(dataframe,2, min)
  
  scaled_data <- scale(dataframe, center = column_mins, scale = column_maxes - column_mins)
  return(data.frame(scaled_data))
  
}





