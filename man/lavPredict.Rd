\name{lavPredict}
\alias{lavPredict}
\alias{lavpredict}
\title{Predict the values of latent or obsered variables.}
\description{
The \code{lavPredict()} function can be used to compute estimated
values for latent variables, and model predicted values for observed
variables.}
\usage{
lavPredict(object, type = "lv", newdata = NULL, method = "EBM",
           se.fit = FALSE, label = TRUE, optim.method = "nlminb")
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{type}{A character string. If \code{"lv"}, estimated values for
the latent variables in the model are computed. If \code{"ov"}, model
predicted values for the observed variables in the model are computed.}
\item{newdata}{An optional data.frame, containing the same variables as
the data.frame used when fitting the model in object.}
\item{method}{A character string. In the linear case (when the indicators are
continuous), the possible options are \code{"regression"} or \code{"Bartlett"}.
In the categorical case, the only option (for now) is \code{"EBM"} for
the Empirical Bayes Modal approach.}
\item{se.fit}{Not used yet.}
\item{label}{Logical. If TRUE, the columns are labeled.}
\item{optim.method}{Character string. Only used in the categorical case. 
If \code{"nlminb"} (the default), the \code{"nlminb()"} function is used 
for the optimization. If \code{"BFGS"}, the \code{"optim()"} function
is used with the BFGS method.}
}
\details{
The \code{predict()} function calls the \code{lavPredict()} function 
with its default options.
}
\seealso{
\code{\link{lavaan}}
}
\examples{
# fit model
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
     
fit <- cfa(HS.model, data=HolzingerSwineford1939)
head(lavPredict(fit))
head(lavPredict(fit, type = "ov"))
}