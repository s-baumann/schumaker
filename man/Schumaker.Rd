% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/Schumaker.R
\name{Schumaker}
\alias{Schumaker}
\title{Create a Schumaker Spline}
\usage{
Schumaker(x, y, ff = "Not-Supplied", Vectorised = TRUE,
  Extrapolation = c("Curve", "Constant", "Linear"))
}
\arguments{
\item{x}{A vector of x coordinates}

\item{y}{A corresponding vector of y coordinates}

\item{ff}{(Optional) A corresponding vector of gradiants at the data points. If not supplied this is estimated.}

\item{Vectorised}{This is a boolean parameter. Set to TRUE if you want to be able to input vectors to the created spline. If you will only input single values set this to FALSE as it is a bit faster.}

\item{Extrapolation}{This determines how the spline function responds when an input is recieved outside the domain of x. The options are "Curve" which outputs the result of the point on the quadratic curve at the nearest interval, "Constant" which outputs the y value at the end of the x domain and "Linear" which extends the spline using the gradiant at the edge of x.}
}
\value{
A list with 3 spline functions. Thee first spline is is for the input points, the second spline is the first derivative of the first spline, the third spline is the second derivative. Each function takes an x value (or vector if Vectorised = TRUE) and outputs the interpolated y value (or relevent derivative).
}
\description{
Create a Schumaker Spline
}
\examples{
x = seq(1,6)
y = log(x)
SSS = Schumaker(x,y, Vectorised = TRUE)
Spline   = SSS[[1]]
SplineD  = SSS[[2]]
Spline2D = SSS[[3]]
xarray = seq(1,6,0.01)
Result = Spline(xarray)
Result2 = SplineD(xarray)
Result3 = Spline2D(xarray)
plot(xarray, Result, ylim=c(-0.5,2))
lines(xarray, Result2, col = 2)
lines(xarray, Result3, col = 3)
}
\references{
Judd (1998). Numerical Methods in Economics. MIT Press
}

