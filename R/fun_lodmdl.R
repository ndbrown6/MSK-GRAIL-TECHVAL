'fun_lodmdl' <- function(df,
			 mdl,
			 grp,
			 ...)
{
	model <- glm(call ~ expected_af, data=df, family = binomial(link = mdl))
	temp.data <- data.frame(expected_af = seq(0.01, max(df$expected_af), 0.01))
	predicted.data <- as.data.frame(predict(model, newdata = temp.data, se = TRUE))
	show.data <- cbind(temp.data, predicted.data) %>%
				 mutate(ymin = model$family$linkinv(fit - 1.96*se.fit),
				        ymax = model$family$linkinv(fit + 1.96*se.fit),
				        yfit = model$family$linkinv(fit),
				        group = grp)
  	return(invisible(show.data))
}
