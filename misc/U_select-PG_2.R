#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

###      F_select: SELECT structural model for paired gear selective analysis 
###                           User routine
###                           Juan Santos - 10.2014
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#input data and starting values:

# example data:

str(x) # this is just and example based on a single haul:

#'data.frame':	106 obs. of  6 variables:
#  $ h   : num  1 1 1 1 1 1 1 1 1 1 ...                   # only one haul included
#$ l   : num  23 25 25.5 26 26.5 27 27.5 28 28.5 29 ...   # fish length with .5 cm precision
#$ ncd1: int  0 0 0 0 0 0 0 0 0 0 ...                     # numbers in codend 1
#$ ncd2: int  1 3 2 2 2 3 2 1 3 4 ...                     # numbers in codend 2
#$ qcd1: num  0.0816 0.0816 0.0816 0.0816 0.0816 ...      # subsampling ration codend 1
#$ qcd2: num  0.049 0.049 0.049 0.049 0.049 ...           # subsampling ratio codend 2


# starting values of the model

str(Start)


#'data.frame':	3 obs. of  4 variables:
#$ sp : num  0.5 0.1 0.9                # starting value for the split parameter, and lower and upper constrains   
#$ l50: int  50 10 90                   # starting value for L50, and lower and upper constrains 
#$ sr : int  10 1 50                    # starting value for SR, and lower and upper constrains 
#$ d  : num  1e+00 1e-04 1e+02          # starting value for D parameter, and lower and upper constrains (only needed for the Richard model)

# if you have experience on non-linear modeling, you will be familiar with chosing tentative starting values. Find the ones which better fits to your specific data



# Modeling the pooled data ------------------------------------------------

# This is only to show you how the model works on pooled data. Actually the results you need to provide needs to account
# for the within- betwen-haul variation vía double bootstrapping.

# The function is called F_select, and you need to enter the following arguments:

        #x  = pooled data

        #lp = range of lengths you want to use for prediction purposes

        # modelo 1= logit
                #2= probit
                #3= gompertz
                #4= Richards
                #0= model selection table: it only outputs fit statistics and AIC values. You should choose the first one (as best model candidate)
        
        #ref_length: optional, used to estimate fishery selectivity indicators based on a reference sizes (e.g. MRCS in European fisheries)
        
        #n.min= minimum number of fish per length class (leave it as one)
        
        # Start= starting values 
        
        # fix_par= tells the model which parameter should be kept fixed during the estimation. fix_par=NULL all parameters will be estimated
                 # it is a good procedure to check if fixing sp=0.5 [fix_par=list("sp",.5)] gives you lower AIC than the same model with fix_par=NULL,
                 # that would indicate equal length-independent catch efficiency among gears. 

        

model_selection<-F_select(x=x,lp=c(10,100),modelo=0,ref_length=NULL,n.min=1,Start=Start,fix_par=NULL)

model_selection_splitFixed<-F_select(x=x,lp=c(10,100),modelo=0,ref_length=NULL,n.min=1,Start=Start,fix_par=list("sp",.5))

rbind(model_selection[1,],model_selection_splitFixed[1,])

   #model rl_fun        sp      l50       sr d   logLik nbetas      aic     aicc
 #      2 probit 0.5201406 53.58236 13.76656 1 1177.971      3 2361.941 2362.176
 #      2 probit 0.5000000 52.52176 13.08257 1 1178.063      2 2360.127 2360.243

#  in the example, the mdoel with fixed sp (second row) provides better aic (and aicc) that the first one, so it should be selected as best candidate


# after deciding which model should be used, simple run it by choosing the proper model (in this case, modelo=2)
best_model_splitFixed<-F_select(x=x,lp=c(10,100),modelo=2,ref_length=NULL,n.min=1,Start=Start,fix_par=list("sp",.5))

#> str(best_model_splitFixed)
#List of 13
#$ fun          : chr "probit"                                # the function you used (modelo 2 = probit)
#$ Betas        : Named num [1:4] 0.5 52.5 13.1 1             # the estimated parameters (here only l50 and sr are estimated, while sp is fixed)
#..- attr(*, "names")= chr [1:4] "split" "l50" "sr" "delta"
#$ nbetas       : int 2                                       # number of betas estimated
#$ l            : num [1:106] 23 25 25.5 26 26.5 27           # vector of  fish lengths observed
#$ L            : num [1:361] 10 10.2 10.5 10.8 11 ...        # prediction grid of lengths
#$ p.phi        : num [1:106] 0 0 0 0 0 0 0 0 0 0 ...         # emprirical catch sharing proportions
#$ selective_ind: NULL                                        # fisheries selectivity indicators. Not used here
#$ phi_check    : num [1:106] 0.00194 0.00377 0 .00442        # predicted catch sharing proportions (using vector of fish lengths observed, and specificvally derived for calculating fit statistics avoiding overinflation caused by subsample ratios)
#$ phi          : num [1:106] 0.00117 0.00227 0.00266         # predicted catch sharing proportions (using vector of fish lengths observed)
#$ PHI          : num [1:361] 5.81e-06 6.54e-06 7.35e-06      # predicted catch sharing proportions (using prediction grid of lengths)
#$ modelhood    : num 1178                                    # model likelihood
#$ aic          : num 2360                                    # aic values
#$ aicc         : num 2360                                    #aicc value


F_check(dat=x,mod=best_model_splitFixed) #gives you residuals, deviance, chi-square p-value (see wileman et al., 1996)

P_select(mod=best_model,nome="",Dir=NULL) # to visualize the data and the fitted curve (can be improved)



# Bootstrap ---------------------------------------------------------------

# Once you have decided (based on AICc) which is the best candidate model to model your data, you have to conduct the boostrap on the original data (stacked haul by haul):

str(X)

#'data.frame':	235 obs. of  6 variables:
#  $ h   : int  2 2 2 2 2 2 2 2 2 5 ...   with h being the haul Id
#$ l   : int  21 24 25 26 28 30 31 34 35 16 ...
#$ ncd1: int  0 0 0 0 0 1 0 0 1 0 ...
#$ ncd2: int  1 1 1 1 1 2 1 1 0 1 ...
#$ qcd1: int  1 1 1 1 1 1 1 1 1 1 ...
#$ qcd2: int  1 1 1 1 1 1 1 1 1 1 ...



#The function to conduct the bootstrap is called B_select_win(). It uses the following arguments:


# X :  data frame with selectivity data by haul, with the structure showed in previous lines:

#Start : same starting values you used in the estimation of the average curve and parameters (using F_select)

# lp : same values you used in the estimation of the average curve and parameters (using F_select)

#n.min :  same values you used in the estimation of the average curve and parameters (using F_select)

#fix_par=list("sp",.5) : same values you used in the estimation of the average curve and parameters (using F_select)

#B : number of pseudosamples in the boostrap scheme 1000 is frequently used

#ncores : The number of cores you want to use

#exportfuns : for windows, you need to tell the machine which functions should be used  when foreach() deliver a iteratioon into a given core. The ffollowing worked for me:

exportfuns<-c(ls()[str_detect(ls(), regex("f_",ignore_case = TRUE))],"funnome","ldply")

# So the bootstrap function should look like that

Bmod<-B_select_win(X=X,Start=Start,modelo=1,lp=c(10,100),n.min=1,fix_par=list("sp",.5),B=10,ncores=4,exportfuns=exportfuns)


#the output of B_select_win is 

# of 10
#$ fun     :   the function you used (as in F_select)
#$ B       :   the number of iterations you used
#$ l       :   vector of  fish lengths observed (as in F_select)
#$ L       :   prediction grid of lengths (as in F_select)
#$ p.phi   :   emprirical catch sharing proportions (as in F_select)
#$ aic     :   aic values of the average model (as in F_select)
#$ aicc    :   aicc values of the average model (as in F_select)
#$ Betas.ci:   Table of parameters with efron confidence intervals 
#$ phi.ci  :   predicted average phi curve  with efron confidence intervals 
#$ betas   :'  1,..,B pseuso-betas estimated during the bootsrap








