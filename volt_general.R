volterra <- function(FreeEl, Core, step = 100, intstart = 0, intend = 10, limits = c(0, 1), method = "rectangle"){
  
  # Possible quadratures for approximating integral
  # "rectangle" is for rule of right rectangles
  # "trapezoid" is for trapezzoid rule
  # "parabol" is for parabol rule (Simpson's rule)
  library("ggplot2")
  
  ms = c("rectangle", "trapezoid", "parabol");
  
  y <- array(0, step + 1)
  h <- (intend - intstart)/step
  dot <- seq(intstart, intend, by = h)
  y[1] = FreeEl(intstart)
  
  if(!(method %in% ms)){stop("Given quadrature argument doesn't fit any of known formulas")}
  
  switch(
    method,
    "rectangle" = {
      # we'll use right rectangle's rule
      for( i in 2:( step + 1 ) )
      {
        sum <- 0
        if( i > 2 )
        {
          for( j in 1:( i - 2 ) )
          {
            sum = sum + Core( dot[i], dot[j+1] ) * y[j + 1] 
          }
        }
        y[i] = ( FreeEl( dot[i] ) + h * sum ) / ( 1 - h * Core( dot[i], dot[i] ) )
      }
    },
    "trapezoid" = {
      for( i in 2:( step + 1 ) )
      {
        sum <- 0
        if( i > 2 )
        {
          for( j in 1:( i - 2 ) )
          {
            sum = sum + Core( dot[i], dot[j + 1]) * y[j + 1] 
          }
        }
        y[i] = ( FreeEl( dot[i] ) + ( h / 2 ) * y[1] * Core( dot[i], intstart) + h * sum ) / ( 1 - ( h / 2 ) * Core( dot[i], dot[i] ) )
      }
    },
    "parabol" = {
      
      y[2] = ( FreeEl( dot[2] ) + ( h / 2 ) * y[1] * Core( dot[2], intstart) ) / ( 1 - ( h / 2 ) * Core( dot[2], dot[2] ) )
      
      for( i in 3:( step + 1 ) )
      {
        sum_even <- 0;
        sum_odd <- 0;
        tryCatch(
          {
            for( j in seq( 1, i - 2, by = 2 ) )
            {
              sum_odd = sum_odd + Core( dot[i], dot[j + 1] ) * y[j + 1]
            }
          }, 
          error = function(e){ sum_odd = 0 }
        )
        tryCatch(
          {
            for(k in seq( 2, i - 2, by = 2 ) )
            {
              sum_even = sum_even + Core( dot[i], dot[k + 1] ) * y[k + 1]
            }
          }, 
          error = function(e){ sum_even = 0 }
        )
        
        if( i %% 2 == 1 )
        {
          y[i] = ( FreeEl( dot[i] ) + ( h / 3 ) * ( Core( dot[i], dot[1] ) * y[1] + 4 * sum_odd + 2 * sum_even ) ) / (1 - ( h / 3 ) * Core( dot[i], dot[i] ) )
        }
        else
        {
          y[i] = ( FreeEl( dot[i] ) + ( h / 2 ) * ( Core( dot[i], dot[1] ) * y[1] + y[2] * Core( dot[i], dot[2] ) ) + ( h / 3 ) * ( Core( dot[i], dot[2] ) * y[2] + 2 * sum_odd + 4 * sum_even ) ) / ( 1 - ( h / 3 ) * Core( dot[i], dot[i] ) )
        }
        
      }
    })
  # Render Graph
  
  if( ( max(y) - min(y) ) < 10 * h ){ limits <- c( ( min(y) - h ) , ( max(y) + h ) ) }
  p <- ggplot( data = data.frame( vals = y, time = dot ), mapping = aes( x = time, y = vals ) )
  p = p + geom_area( fill = "lightblue", color = "darkblue" )
  p = p + labs( x = "Time", y = "", title = "Walther Equation" )
  p = p + theme( plot.title = element_text( hjust = 0.5, face = "bold" ) )
  p = p + ylim( limits )
  print(p)
  res = list(values = y, plot = p)
  return( res )
}
