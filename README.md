# Project to estimate Chinook salmon vital rates
This repository contains an example model for estimating Chinook salmon vital rates and a simulation exercise to test model identifiability, convergence, and performance. This model, MLB 4, assumes additive age and cohort/year effects on the complementary log-log scale for exploitation, maturation, and natural mortality rates. 

To start open 'mlb4.r' and input the type of data that should be fit (stochastic or deterministic/noise free), the number of data sets (n.data.sets) to fit (can be greater than 1 if fitting stochastic data), and number of times to refit each data set using different starting points (n.iter). To plot all solutions of a single data set (randomly selected if more than 1), run 'MLB4_plot.1dataset.gen.r'. This script currently plots all solutions but can easily be adjusted to only plot solutions that satisfy a given convergence criteria. 