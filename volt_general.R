volt <- function(FreeEl, Core, step = 100, intstart = 0, intend = 10, limits = c(0, 1), method = "rectangle", method_global = "sum", eps = 0.1){
  
  # Possible quadratures for approximating integral
  # "rectangle" is for rule of right rectangles
  # "trapezoid" is for trapezzoid rule
  # "parabol" is for parabol rule (Simpson's rule)
  # "global_method" defines global rule for numeric solution - approximation method or method of finite sums
  if ( !require( "ggplot2" ) )
  {
    install.packages( "ggplot2" )
    if( !require( "ggplot2") ) stop("Package not found")
  }
  if ( !require( "Deriv" ) )
  {
    install.packages( "Deriv" )
    if( !require( "Deriv") ) stop("Package not found")
  }
  ms <- c("rectangle", "trapezoid", "parabol");
  
  y <- array(0, step + 1)
  h <- (intend - intstart)/step
  dot <- seq(intstart, intend, by = h)
  z <- array(0, step + 1)
  z[1] =  y[1] = FreeEl(intstart)
  

  switch(method_global,
    "sum" = {
        f1 <- Deriv(FreeEl, x = "g")
        f2 <- Deriv(f1, x = "g")
        R <- array(0, step + 1)
        err <- array(0, step + 1)
        
        if(!(method %in% ms)){stop("Given quadrature argument doesn't fit any of known formulas")}
        
        switch(
          method,
          "rectangle" = {
            # we'll use right rectangle's rule
            for( i in 2:( step + 1 ) )
            {
              R[i] = ( ( h ** 2 ) / 24 ) * (dot[i] - dot[i - 1]) * optimise(f2, interval = c(dot[i - 1], dot[i]), maximum = TRUE )$objective
              sum <- 0
              sum_e <- 0
              if( i > 2 )
              {
                for( j in 1:( i - 2 ) )
                {
                  sum = sum + Core( dot[i], dot[j + 1] ) * y[j + 1] 
                  sum_e = sum_e + Core( dot[i], dot[j + 1]) * err[j + 1]
                }
              }
              y[i] = ( FreeEl( dot[i] ) + h * sum ) / ( 1 - h * Core( dot[i], dot[i] ) )
              err[i] = ( ( sum_e + R[i] ) / ( 1 - Core( dot[i], dot[i] ) ) )
            }
          },
          "trapezoid" = {
            for( i in 2:( step + 1 ) )
            {
              R[i] = ( ( h ** 2 ) / 12 ) * (dot[i] - dot[i - 1]) * optimise(f2, interval = c(dot[i - 1], dot[i]), maximum = TRUE )$objective
              sum <- 0
              sum_e <- 0
              if( i > 2 )
              {
                for( j in 1:( i - 2 ) )
                {
                  sum = sum + Core( dot[i], dot[j + 1]) * y[j + 1] 
                  sum_e = sum_e + Core( dot[i], dot[j + 1]) * err[j + 1]
                }
              }
              y[i] = ( FreeEl( dot[i] ) + ( h / 2 ) * y[1] * Core( dot[i], intstart) + h * sum ) / ( 1 - ( h / 2 ) * Core( dot[i], dot[i] ) )
              err[i] = ( ( sum_e + R[i] ) / ( 1 - Core( dot[i], dot[i] ) ) )
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
        },
        "approx" = {
          s <- 1
          for( k in 2:( step + 1 ) )
          {
            z[k] = FreeEl( ( k - 1 ) * h )
          }
          while( abs( s ) > eps )
          {
            s <- 0
            for( i in 2:( step + 1 ) )
            {
              sum <- 0
              if( i > 2 )
              {
                for( j in 1:(i - 2) )
                {
                  sum = sum + Core( ( i - 1 ) * h, j * h) * z[j + 1]
                }
              }
              y[i] = ( FreeEl( ( i - 1 ) * h ) + ( h / 2 ) * z[1] * Core( (i - 1) * h, intstart ) + h * sum + ( h / 2 ) * Core( ( i - 1 ) * h, ( i - 1 ) * h ) * z[i] )
              s = s + abs( y[i] - z[i] )
              z[i] = y[i]
            }
          }
          err = eps
          })
  # Render Graph
  
  if( ( max(y) - min(y) ) < 10 * h ){ limits <- c( ( min(y) - h ) , ( max(y) + h ) ) }
  p <- ggplot( data = data.frame( vals = y, time = dot ), mapping = aes( x = time, y = vals ) )
  p = p + geom_area( fill = "lightblue", color = "darkblue" )
  p = p + labs( x = "Time", y = "", title = "Volterra Equation" )
  p = p + theme( plot.title = element_text( hjust = 0.5, face = "bold" ) )
  p = p + ylim( limits )
  print(p)
  res = list(values = y, plot = p, error = max(err) )
  return( res )
}
