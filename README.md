# R-code-for-signal-lasso-type
This repository contains the R code for signal-lasso-type methods
main_slt.cpp is an exporting C++ function to R. You can source this function into an R session using the 
Rcpp::sourceCpp function (or via the Source button on the editor toolbar). Learn more about Rcpp at:
http://www.rcpp.org/
http://adv-r.had.co.nz/Rcpp.html
http://gallery.rcpp.org/

main_slt.R is main function for calculation of estimator of signal-lasso-type method, which including signal lasso, adaptive signal lasso, signal lasso with non-convex penalty (product and minimum penalty) 

simulation_slt gives the R-code or functions for calculating various shringkage methods,including lasso, adaptive lasso, 
SCAD, MCP, ElasticNet, signal lasso, adaptive signal lasso, signal lasso with non-convex penalty for purpose of comparisons.

simu_exam.R gives the example of calculation in linear regression model. 
