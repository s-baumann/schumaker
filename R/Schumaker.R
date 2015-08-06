
#' Create a Schumaker Spline
#' @export
#' @param tt A vector of x coordinates
#' @param FF A corresponding vector of y coordinates
#' @param ff (Optional) A corresponding vector of gradiants at the data points. If not supplied this is estimated.
#' @param Vectorised This is a boolean parameter. Set to TRUE if you want to be able to input vectors to the created spline. If you will only input single values set this to FALSE as it is a bit faster.
#' @param Extrapolation This determines how the spline function responds when an input is recieved outside the domain of tt. The options are "Curve" which outputs the result of the point on the quadratic curve at the nearest interval, "Constant" which outputs the FF value at the end of the tt domain and "Linear" which extends the spline using the gradiant at the edge of tt.
#'
#' @return A list with 3 spline functions. Thee first spline is is for the input points, the second spline is the first derivative of the first spline, the third spline is the second derivative. Each function takes an x value (or vector if Vectorised = TRUE) and outputs the interpolated y value (or relevent derivative).
#' @examples
#' tt = c(0,1,2,3,4,5,6)
#' FF = c(2,1, 0.59, 0.27, 0.25, -0.23, -0.45)
#'
#' SSS = Schumaker(tt,FF, Vectorised = TRUE)
#' Spline   = SSS[[1]]
#' SplineD  = SSS[[2]]
#' Spline2D = SSS[[3]]
#'
#' xarray = seq(0,6,0.1)
#' Res = Spline(xarray)
#' Res2 = SplineD(xarray)
#' Res3 = Spline2D(xarray)
#'
#' plot(xarray, Res, ylim=c(-2,5))
#' lines(xarray, Res2, col = 2)
#' lines(xarray, Res3, col = 3)

Schumaker <- function(tt,FF, ff = "Not-Supplied", Vectorised = FALSE, Extrapolation = c("Curve", "Constant", "Linear")){
  Extrapolation = Extrapolation[1]
  if (!(Extrapolation %in% c("Constant", "Linear", "Curve"))){stop("The extrapolation parameter defines what the function returns when evaluated
                                                                      outside the domain of the interpolation data. \n Choose 'Constant' for constant
                                                                      extrapolation. This returns the value at the nearest edge of the domain. \n 'Linear'
                                                                      extends out a line from the edge of the domain with a slope of the gradiant at
                                                                      that point. \n 'Curve' extrapolation uses the parabolic arc for the last interval.")}

# Schumaker shape-preserving quadratic interpolation spline.

n = length(tt)

if (ff == "Not-Supplied"){
  # Judd (1998), page 233, second last equation
  L = sqrt( (tt[2:n]-tt[1:(n-1)])^2 + (FF[2:n]-FF[1:(n-1)])^2)
  # Judd (1998), page 233, last equation
  d = (FF[2:n]-FF[1:(n-1)])/(tt[2:n]-tt[1:(n-1)])
  # Judd (1998), page 234, Eqn 6.11.6
  Conditionsi = (d[1:(n-2)]*d[2:(n-1)] > 0)
  MiddleSiwithoutApplyingCondition = (L[1:(n-2)]*d[1:(n-2)]+L[2:(n-1)] * d[2:(n-1)]) / (L[1:(n-2)]+L[2:(n-1)])
  sb = Conditionsi * MiddleSiwithoutApplyingCondition
  # Judd (1998), page 234, Second Equation line plus 6.11.6 gives this array of slopes.
  ff = c(((-sb[1]+3*d[1])/2),  sb,  ((3*d[n-1]-sb[n-2])/2))
}

NumberOfIntervalsWithKnots = 2*(n-1)

Intervals = 1:(n-1)
IntervalTab = data.frame(IntervalNum = sort(rep(Intervals,2)),
                                    SubIntervalNum = rep(c(1,2),n-1),
                                    StartOfInterval = numeric(NumberOfIntervalsWithKnots),
                                    EndOfInterval = numeric(NumberOfIntervalsWithKnots)
                                    )

Evals = do.call(rbind, lapply(Intervals, function(IntervalNum) SchumakerIndInterval(c(FF[IntervalNum], FF[IntervalNum+1]), c(ff[IntervalNum], ff[IntervalNum+1]), c(tt[IntervalNum], tt[IntervalNum+1]))))

IntervalTab = cbind(IntervalTab, Evals)
rm(Evals)


IntervalTab[IntervalTab$SubIntervalNum == 1, "StartOfInterval"] = tt[1:n-1]
IntervalTab[IntervalTab$SubIntervalNum == 2, "EndOfInterval"] = tt[2:n]
IntervalTab[IntervalTab$SubIntervalNum == 2, "StartOfInterval"] = IntervalTab$tsi[IntervalTab$SubIntervalNum == 2]
IntervalTab[IntervalTab$SubIntervalNum == 1, "EndOfInterval"] = IntervalTab$tsi[IntervalTab$SubIntervalNum == 1]
# This gets rid of empty intervals. The 1e-10 is there in case of numerical imprecision.
IntervalTab <<- IntervalTab[which(IntervalTab$EndOfInterval + 1e-10 > IntervalTab$StartOfInterval),]

if ((Extrapolation %in% c("Constant", "Linear"))){
  Botx = min(IntervalTab$StartOfInterval)
  Boty   = FF[1]

  BotInt = findInterval(Botx, IntervalTab$StartOfInterval)
  BotB = IntervalTab[BotInt ,c("B")]
  BotC   = Boty - BotB

  cat("Test Test, \n")

    if (Extrapolation == "Constant"){

      BotB = 0
      BotC   = Boty
    }

  BotRow = data.frame(IntervalNum = 0, SubIntervalNum = 0, StartOfInterval = Botx-1, EndOfInterval = Botx, tsi = 0, C = 0, B =  BotB, A = BotC)


  Topx = max(IntervalTab$EndOfInterval)
  Topy   = FF[n]

  TopInt = findInterval(Topx, IntervalTab$StartOfInterval)

  TopB = IntervalTab[TopInt ,c("B")]
  TopC   = Topy

  if (Extrapolation == "Constant"){
    TopB = 0
    TopC   = Topy
  }

  TopRow = data.frame(IntervalNum = 0, SubIntervalNum = 0,StartOfInterval = Topx, EndOfInterval = Topx + 1, tsi = 0, C = 0, B =  TopB, A = TopC)

  IntervalTab = rbind(BotRow, IntervalTab, TopRow)
}

# It is important use individual vectors and matrices rather than datatables or data.frames for speed
IntStarts = c(IntervalTab$StartOfInterval, Inf)
SpCoefs = data.matrix(IntervalTab[,c("C", "B", "A")])

# This is the end spline which looks up the correct interval and evaluates with the correpsonding coefficients
Spline0 = ppmak(IntStarts, SpCoefs, Vectorised )
Spline1 = ppmakDeriv(IntStarts, SpCoefs, Vectorised )
Spline2 = ppmak2Deriv(IntStarts, SpCoefs, Vectorised )
# This just boosts the speed of evaluation by 5ish percent. Not essential.
CompiledSpline0 = compiler::cmpfun(Spline0)
CompiledSpline1 = compiler::cmpfun(Spline1)
CompiledSpline2 = compiler::cmpfun(Spline2)
return(list(Spline = CompiledSpline0, DerivativeSpline = CompiledSpline1, SecondDerivativeSpline = CompiledSpline2 ))
}
