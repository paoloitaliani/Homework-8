# Homework-8
Weekly assignment about dimension reduction for linear regression

## Exercise 1

On the Virtuale page of the course (under Data Sets) you find the dataset
LBPnumerical.dat. The 773 observations are lower back pain patients. The Variable y is the Roland-Morris
score, a measurement of the severity of back pain after two weeks of treatment on a scale
between 0 and 100. It is of interest to predict this from the 102 x-variables, which all
contain information about the patients before their treatment was started. This dataset
is a subset (with some observations and variables removed because of too many missing
values, and some imputed values for further missing values, so that no missing values are
left in the data) of the dataset documented on

http://ifcs.boku.ac.at/repository/challenge2/

so you can go there if you like to find out more about the variables, which is not strictly
necessary for this exercise.
The main interest here is to construct a good prediction model for the Roland-Morris
score. Find such a model from a comparison of approaches that should involve at least the
lasso, principal component regression, and a comparison by double cross-validation. Also
demonstrate that (or whether) your model is better than a model that predicts y from the
mean only. You can try out more if you want.


## Exercise 2 

For the Covid-data, try to find a model that is substantially better than prediction by mean only using the methods that you have learnt up to now, using a logtransformed y-variable (note that for comparison you also need to assess the mean-only
model using the transformed y). If you think it is useful, you can also transform x-variables.
