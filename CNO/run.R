rm(list = ls())
library(CellNOptR)
library(Rgraphviz)
tic = Sys.time()

# Paper: 
# https://www.dropbox.com/s/94v6wouuszsycg4/Terfve%20et%20al.%20-%202012%20-%20CellNOptR%20a%20flexible%20toolkit%20to%20train%20protein%20sig.pdf?dl=0

# 실험 데이터를 불러들이고, 불러들인 데이터가 어떤 모양인지 시각화합니다. 
ToyDataRail = CNOlist("ToyDataMMB.csv")
pdf("datarail.pdf")
plotCNOlist(ToyDataRail)
dev.off()

# 네트워크 구조를 불러들입니다. 이것은 실험데이터, 문헌데이터를 이용하여 수립한 네트워크일 것입니다. 
Toy<-readSIF("ToyModelMMB2.sif")
pdf("Toy.pdf")
plotModel(Toy, ToyDataRail)
dev.off()

# 2. CUT NCNO COMPONENTS FROM THE MODEL 
# Finds the indices, in the model fields, of the species that 
# are measured/inhibited/stimulated. This function finds the indices, 
# in the model fields, of the species that are measured/inhibited/
# stimulated. It looks for their position in model$namesSpecies which has 
# the same order as the rows of interMat and notMat, and therefore these 
# indexes can be used there as well.
iToy<-indexFinder(ToyDataRail, Toy, verbose=TRUE)

# Find the indexes of the non-observable and non controllable species
iToyNCNO<-findNONC(Toy, iToy, verbose=TRUE)

# Cuts the non-observable/non-controllable species from the model
ToyNCNOcut<-cutNONC(Toy, iToyNCNO)

iToyNCNOcut<-indexFinder(ToyDataRail, ToyNCNOcut, verbose=TRUE)

pdf("Toy-NCNOcut.pdf")
plotModel(ToyNCNOcut, ToyDataRail)
dev.off()

# 3. COMPRESSING THE MODEL 
ToyNCNOcutComp<-compressModel(ToyNCNOcut, iToyNCNOcut)
iToyCNCOcutCOmp<-indexFinder(ToyDataRail, ToyNCNOcutComp) 

pdf("Toy-NCNOcutComp.pdf")
plotModel(ToyNCNOcutComp, ToyDataRail)
dev.off()

# 4. EXPANDING THE GATES 
ToyNCNOcutCompGatexp<-expandGates(ToyNCNOcutComp, maxInputsPerGate=3)
pdf("Toy-NCNOcutCompGatexp.pdf")
plotModel(ToyNCNOcutCompGatexp, ToyDataRail)
dev.off()

# 5. OPTIMIZING THE MODEL

resToy<-residualError(ToyDataRail) 
initBstring<-rep(1, length(ToyNCNOcutCompGatexp$reacID))

ToyT1opt<-gaBinaryT1(CNOlist=ToyDataRail, model=ToyNCNOcutCompGatexp, 
    initBstring=initBstring, verbose=TRUE
    )

pdf("GA-result.pdf")
cutAndPlot(model=ToyNCNOcutCompGatexp, 
    bStrings=list(ToyT1opt$bString), CNOlist=ToyDataRail
    )
dev.off()

pdf("GA-fitness.pdf")
plotFit(optRes=ToyT1opt)
dev.off()

# 6. DISPLAY THE OPTIMIZED MODEL NETWORK 
pdf("Toy-NCNOcutCompGatexp-Opt.pdf")
plotModel(ToyNCNOcutCompGatexp,ToyDataRail,bString=ToyT1opt$bString)
dev.off()

bs = mapBack(ToyNCNOcutCompGatexp, Toy, ToyT1opt$bString)
pdf("Toy-NCNOcutCompGatexp-Opt-Mapback.pdf")
plotModel(Toy,ToyDataRail,bs,compressed=Toy$speciesCompressed)
dev.off()

toc <- Sys.time() - tic 
print(toc)

# 7. write 
writeScaffold(
    modelOriginal=Toy,
    modelComprExpanded=ToyNCNOcutCompGatexp,
    optimResT1=ToyT1opt,
    optimResT2=NA,
    CNOlist=ToyDataRail
    )

dotObj<-agread("Scaffold.dot",layoutType="dot",layout=TRUE)
toFile(dotObj,layoutType="dot",filename="Scaffold.ps",fileType="ps")

writeNetwork(
    modelOriginal=Toy,
    modelComprExpanded=ToyNCNOcutCompGatexp,
    optimResT1=ToyT1opt,
    optimResT2=NA,
    CNOlist=ToyDataRail
    )

dotObj<-agread("PKN.dot",layoutType="dot",layout=TRUE)
toFile(dotObj,layoutType="dot",filename="PKN.ps",fileType="ps")
