# CNOR시뮬레이션 테스트이다. 
#source("http://bioconductor.org/biocLite.R")
#biocLite("CNORode")
#biocLite("MEIGOR")
#install.packages("MEIGOR_0.99.3_svn2444.tar.gz",type="source")

library(CNORode)
# model=readSIF(system.file("doc", "ToyModelMMB_FeedbackAnd.sif", 
#     package="CNORode"));
# cno_data=readMIDAS(system.file("doc", "ToyModelMMB_FeedbackAnd.csv", 
#     package="CNORode"));
model=readSIF("ToyModelMMB_FeedbackAnd.sif");
cno_data=readMIDAS("ToyModelMMB_FeedbackAnd.csv");
cnolist=makeCNOlist(cno_data,subfield=FALSE);
ode_parameters=createLBodeContPars(
    model,
    LB_n = 1,
    LB_k = 0.1,
    LB_tau = 0.01,
    UB_n = 5,
    UB_k = 0.9,
    UB_tau = 10, 
    default_n = 3,
    default_k = 0.5,
    default_tau = 1,
    opt_n = TRUE,
    opt_k = TRUE,
    opt_tau = TRUE,
    random = FALSE
    )


print(ode_parameters)

modelSim=plotLBodeModelSim(
    cnolist,
    model,
    ode_parameters,
    timeSignals=seq(0,2,0.5)
    );


initial_pars=createLBodeContPars(
    model,
    LB_n = 1, 
    LB_k = 0.1,
    LB_tau = 0.01, 
    UB_n = 5, 
    UB_k = 0.9, 
    UB_tau = 10, 
    random = TRUE
    )
#Visualize initial solution
simulatedData=plotLBodeFitness(cnolist, model,initial_pars)
paramsGA = defaultParametersGA()
paramsGA$maxStepSize = 1
paramsGA$popSize = 50
paramsGA$iter = 100
paramsGA$transfer_function = 2
opt_pars=parEstimationLBode(
    cnolist,
    model,
    ode_parameters=initial_pars,
    paramsGA=paramsGA
    )
#Visualize fitted solution
simulatedData=plotLBodeFitness(cnolist, model,ode_parameters=opt_pars)

simulatedData=plotLBodeFitness(cnolist, model,initial_pars)

simulatedData=plotLBodeFitness(cnolist, model,ode_parameters=opt_pars)

library(MEIGOR)

f_hepato<-getLBodeContObjFunction(cnolist, model, initial_pars, indices=NULL,
    time = 1, verbose = 0, transfer_function = 2, reltol = 1e-05, atol = 1e-03,
    maxStepSize = Inf, maxNumSteps = 1e4, maxErrTestsFails = 50, nan_fac = 1)

n_pars=length(initial_pars$LB);

problem<-list(f=f_hepato, x_L=initial_pars$LB[initial_pars$index_opt_pars],
    x_U=initial_pars$UB[initial_pars$index_opt_pars]);

#Source a function containing the options used in the CeSSR publication
source(system.file("benchmarks","get_paper_settings.R",package="MEIGOR"))
#Set max time as 20 seconds per iteration
opts<-get_paper_settings(20);

Results<-CeSSR(problem,opts,Inf,Inf,3,TRUE,
    global_save_list=c('cnolist','model', 'initial_pars'))

