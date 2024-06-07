# S3 pagfl

    Call:
    pagfl(formula = y ~ a + b, data = data, n_periods = 150, lambda = 5)
    
    Balanced panel: N = 20, T = 150, obs = 3000
    
    Convergence reached:
    TRUE (542 iterations)
    
    Information criterion:
          IC   lambda 
    1.045192 5.000000 
    
    Residuals:
         Min       1Q   Median       3Q      Max 
    -3.10271 -0.67485 -0.01057  0.68122  3.23017 
    
    3 groups:
     1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 
     1  1  1  2  1  1  2  2  2  3  3  2  3  1  3  3  3  1  2  2 
    
    Coefficients:
                    a        b
    Group 1 -0.98847  1.54126
    Group 2 -5.05014 -1.02301
    Group 3  0.22479  1.54219
    
    Residual standard error: 0.99552 on 2978 degrees of freedom
    Mean squared error 0.9838
    Multiple R-squared: 0.91956, Adjusted R-squared: 0.91899 

---

    Code
      estim
    Output
      Groups: 3 
      
      Call:
      pagfl(formula = y ~ a + b, data = data, n_periods = 150, lambda = 5)
      
      Coefficients:
                     a        b
      Group 1 -0.98847  1.54126
      Group 2 -5.05014 -1.02301
      Group 3  0.22479  1.54219

# S3 tv_pagfl

    Call:
    tv_pagfl(formula = y ~ X1, data = data, n_periods = 100, lambda = 7)
    
    Balanced panel: N = 10, T = 100, obs = 1000
    
    Convergence reached:
    TRUE (1546 iterations)
    
    Information criterion:
         IC  lambda 
    1.79141 7.00000 
    
    Residuals:
         Min       1Q   Median       3Q      Max 
    -3.16362 -0.65668 -0.00871  0.72063  2.79348 
    
    5 groups:
     1  2  3  4  5  6  7  8  9 10 
     1  1  2  3  4  5  5  3  5  3 
    
    Residual standard error: 1.02363 on 980 degrees of freedom
    Mean squared error 1.02686
    Multiple R-squared: 0.81416, Adjusted R-squared: 0.81055 

---

    Code
      estim
    Output
      Groups: 5 
      
      Call:
      tv_pagfl(formula = y ~ X1, data = data, n_periods = 100, lambda = 7)

# S3 tv_pagfl const coef unbalanced

    Call:
    tv_pagfl(formula = y ~ X + a, data = df, index = c("i_index", 
        "t_index"), lambda = 25, const_coef = "a")
    
    Unbalanced panel: N = 10, T = 60-77, obs = 698
    
    Convergence reached:
    TRUE (5497 iterations)
    
    Information criterion:
          IC   lambda 
     1.48316 25.00000 
    
    Residuals:
         Min       1Q   Median       3Q      Max 
    -3.57228 -0.66901  0.04044  0.67613  3.03945 
    
    2 groups:
     1  2  3  4  5  6  7  8  9 10 
     1  1  2  1  1  2  2  1  2  1 
    
    Constant coefficients:
                   a
    Group 1 0.01486
    Group 2 0.01400
    
    Residual standard error: 1.08735 on 677 degrees of freedom
    Mean squared error 1.14676
    Multiple R-squared: 0.7833, Adjusted R-squared: 0.7769 

---

    Code
      estim
    Output
      Groups: 2 
      
      Call:
      tv_pagfl(formula = y ~ X + a, data = df, index = c("i_index", 
          "t_index"), lambda = 25, const_coef = "a")

