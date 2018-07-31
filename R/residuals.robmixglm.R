residuals.robmixglm <-
  function(object,
           type = c("deviance", "pearson"),
           ...)
  {
    deviance.residuals <- function(){
      calc.binomial <- function(){
        n <- object$model$Y[,1]+object$model$Y[,2]
        sign(object$model$Y[,1]-n*fitted(object))*sqrt(
          2*(object$model$Y[,1]*log(object$model$Y[,1]/(n*fitted(object)))+
             object$model$Y[,2]*log(object$model$Y[,2]/(n*(1-fitted(object))))))
      }
      res <- switch(object$family,
            gaussian=(object$model$Y-fitted(object))/sqrt(coef(object$fit)[length(coef(object$fit))]),
            binomial=calc.binomial(),
            poisson=sign(object$model$Y-fitted(object))*sqrt(2*(object$model$Y*log(object$model$Y/fitted(object))-(object$model$Y-fitted(object)))),
            gamma=sign(object$model$Y-fitted(object))*sqrt(2*(-log(object$model$Y/fitted(object))+(object$model$Y-fitted(object))/fitted(object))),
            truncpoisson=sign(object$model$Y-fitted(object))*sqrt(2*(object$model$Y*log(object$model$Y/fitted(object))-(object$model$Y-fitted(object))))
      )
      res
    }
    pearson.residuals <- function(){
      res <- switch(object$family,
           gaussian=(object$model$Y-fitted(object))/sqrt(coef(object$fit)[length(coef(object$fit))]),
           binomial=(object$model$Y[,1]-(object$model$Y[,1]+object$model$Y[,2])*fitted(object))/
             sqrt((object$model$Y[,1]+object$model$Y[,2])*fitted(object)*(1-fitted(object))),
           poisson=(object$model$Y-fitted(object))/sqrt(fitted(object)),
           gamma=(object$model$Y-fitted(object))/fitted(object),
           truncpoisson=sign(object$model$Y-fitted(object))/(sqrt(fitted(object)*exp(fitted(object))*(-1+exp(fitted(object))-fitted(object)))/(exp(fitted(object))-1))
      )
      res
    }
    type <- match.arg(type)
     res <- switch(type,
           deviance=deviance.residuals(),
           pearson=pearson.residuals()
     )
    res
  }