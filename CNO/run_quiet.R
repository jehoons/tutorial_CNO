rm(list = ls())
library(CellNOptR)
library(Rgraphviz)
tic = Sys.time()

ToyDataRail = CNOlist("ToyDataMMB.csv")

Toy<-readSIF("ToyModelMMB2.sif")

iToy<-indexFinder(ToyDataRail, Toy, verbose=TRUE)

iToyNCNO<-findNONC(Toy, iToy, verbose=TRUE)

ToyNCNOcut<-cutNONC(Toy, iToyNCNO)

iToyNCNOcut<-indexFinder(ToyDataRail, ToyNCNOcut, verbose=TRUE)

ToyNCNOcutComp<-compressModel(ToyNCNOcut, iToyNCNOcut)

iToyCNCOcutCOmp<-indexFinder(ToyDataRail, ToyNCNOcutComp) 

ToyNCNOcutCompGatexp<-expandGates(ToyNCNOcutComp, maxInputsPerGate=3)

resToy<-residualError(ToyDataRail) 

initBstring<-rep(1, length(ToyNCNOcutCompGatexp$reacID))

ToyT1opt<-gaBinaryT1(CNOlist=ToyDataRail,model=ToyNCNOcutCompGatexp, 
    initBstring=initBstring,verbose=FALSE)

# 7. write 
writeScaffold(modelOriginal=Toy,modelComprExpanded=ToyNCNOcutCompGatexp,
    optimResT1=ToyT1opt,optimResT2=NA,CNOlist=ToyDataRail)

