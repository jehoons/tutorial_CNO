# CellNetOpt

`CellNetOpt`는 불리언모델의 최적화를 위한 툴이다
[(paper)](https://www.dropbox.com/s/94v6wouuszsycg4/Terfve%20et%20al.%20-%202012%20-%20CellNOptR%20a%20flexible%20toolkit%20to%20train%20protein%20sig.pdf?dl=0),[(manual)](https://www.bioconductor.org/packages/devel/bioc/manuals/CellNOptR/man/CellNOptR.pdf).

최적화 결과는 `hypergraph`형식을 가지는데, 예를 들면, 노드 이름이 `and1..n`을 포함한다는 것을 의미한다. 이들을 이용하여 로직방정식을 만들 수 있다. 어떻게? 먼저 and의 소스노드들은 A and B와 같은 형식으로 로직함수를 만들어 준다. . 이후 로직함수로 축약해 준다.

예를 들어서 아래와 같은 형식으로 하이퍼그래프형식의 출력이 주어지는 경우,
```
A 1 and1
B 1 and1
and1 1 C
D 1 C
```

1단계
```
and1 = A and B
C = and1 or D
```

2단계
```
C = (A and B) or D
```

## Install

**CellNOptR 설치하기**

```R
# 컬러풀한 R의 command-line을 원한다면 colorout 을 설치하면 된다.
git clone https://github.com/jalvesaq/colorout.git
R CMD INSTALL colorout

# RCurl을 설치하는 중 에러가 나는 경우
http://askubuntu.com/questions/359267/cannot-find-curl-config-in-ubuntu-13-04

# try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("CellNOptR")
```

**도움말 보기**

```R
browseVignettes("CellNOptR")
```

**파이썬 레퍼 설치**

```
pip install cellnopt.wrapper
```

파이썬 레퍼는 내부적으로 rpy2를 실행합니다. rpy2 를 내부적으로 실행하는 과정에서 레포지토리를 찾지 못한다는 이야기가 나오면, `lib/python2.7/site-packages/rtools/package.py` 코드의  `cran_repos`변수를 "https://cloud.r-project.org/" 값으로 설정해 준다.

## Usage

다음 명령들을 순차적으로 실행해 보자.
```R
rm(list = ls())
library(CellNOptR)
library(Rgraphviz)
tic = Sys.time()
```

먼저 [experimental data](https://github.com/jehoons/lab/blob/master/codebase_2/cellnetopt/ToyDataMMB.csv)를 불러들인다. 데이터의 형식은 [DataRail](https://www.ncbi.nlm.nih.gov/pubmed/18218655)이라는 형식을 따라서 기록된다. 그리고 불러들인 데이터가 어떤 모양인지 시각화한다.
```R
ToyDataRail = CNOlist("ToyDataMMB.csv")
pdf("datarail.pdf")
plotCNOlist(ToyDataRail)
dev.off()
```

`sif`형식으로 정의된 [네트워크](https://github.com/jehoons/lab/blob/master/codebase_2/cellnetopt/ToyModelMMB2.sif)를 불러들인다. 이것은 실험데이터, 문헌데이터를 이용하여 수립한 네트워크일 것이다. 그리고 이것을 시각화한다.
```R
Toy<-readSIF("ToyModelMMB2.sif")
pdf("Toy.pdf")
plotModel(Toy, ToyDataRail)
dev.off()
```

`NCNO COMPONENTS`는 관찰할수도 없고 제어할수도 없는 요소를 의미한다. 이것들을 찾아서 제거한다.

Finds the indices, in the model fields, of the species that are measured/inhibited/stimulated. This function finds the indices, in the model fields, of the species that are measured/inhibited/ stimulated. It looks for their position in model$namesSpecies which has the same order as the rows of interMat and notMat, and therefore these indexes can be used there as well.

그러기 위해서는 우선 모델 필드에서 측정되고, 저해되고, 자극된 분자종에 대한 인덱스를 찾는다.
```R
iToy<-indexFinder(ToyDataRail, Toy, verbose=TRUE)
```

다음으로는 NONC를 찾는다. Find the indexes of the non-observable and non controllable species
```R
iToyNCNO<-findNONC(Toy, iToy, verbose=TRUE)
```

그리고 NONC를 제거한 모델을 얻는다.
```R
ToyNCNOcut<-cutNONC(Toy, iToyNCNO)
iToyNCNOcut<-indexFinder(ToyDataRail, ToyNCNOcut, verbose=TRUE)
pdf("Toy-NCNOcut.pdf")
plotModel(ToyNCNOcut, ToyDataRail)
dev.off()
```

다음으로는 모델을 압축한다.
```R
ToyNCNOcutComp<-compressModel(ToyNCNOcut, iToyNCNOcut)
iToyCNCOcutCOmp<-indexFinder(ToyDataRail, ToyNCNOcutComp)
pdf("Toy-NCNOcutComp.pdf")
plotModel(ToyNCNOcutComp, ToyDataRail)
dev.off()
```

우리는 이전에 `sif` 네트워크 형식을 불렀다. sif 네트워트로부터 게이트함수로 확장한다.
```R
ToyNCNOcutCompGatexp<-expandGates(ToyNCNOcutComp, maxInputsPerGate=3)
pdf("Toy-NCNOcutCompGatexp.pdf")
plotModel(ToyNCNOcutCompGatexp, ToyDataRail)
dev.off()
```

GA를 이용하여 모델을 최적화한다. 이를 위해서 우선 초기화하고, 그다음 GA를 실행한다.
```R
resToy<-residualError(ToyDataRail)
initBstring<-rep(1, length(ToyNCNOcutCompGatexp$reacID))
ToyT1opt<-gaBinaryT1(CNOlist=ToyDataRail, model=ToyNCNOcutCompGatexp,
    initBstring=initBstring, verbose=TRUE
    )
```

그리고 트레이닝 결과를 시각화한다.
```R
pdf("GA-result.pdf")
cutAndPlot(model=ToyNCNOcutCompGatexp,
    bStrings=list(ToyT1opt$bString), CNOlist=ToyDataRail
    )
dev.off()
pdf("GA-fitness.pdf")
plotFit(optRes=ToyT1opt)
dev.off()
```

최적화된 모델의 구조를 시각화한다.
```R
pdf("Toy-NCNOcutCompGatexp-Opt.pdf")
plotModel(ToyNCNOcutCompGatexp,ToyDataRail,bString=ToyT1opt$bString)
dev.off()
bs = mapBack(ToyNCNOcutCompGatexp, Toy, ToyT1opt$bString)
pdf("Toy-NCNOcutCompGatexp-Opt-Mapback.pdf")
plotModel(Toy,ToyDataRail,bs,compressed=Toy$speciesCompressed)
dev.off()
toc <- Sys.time() - tic
print(toc)
```

결과를 기록한다. [writeScaffold](https://rdrr.io/bioc/CellNOptR/man/writeScaffold.html)함수를 이용.
```R
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
```

## Idea - 구조를 변경하면서 최적화 하기.
일반적으로 초기의 모델 구조가 실험데이터에 잘 부합하는 경우는 잘 없을 것이다. 이를 위해서, 이 코드를 실행하기 전에 각 분자들의 PPI 상호작용에 대한 정보(예: STRING)를 이용하여 가능한 네트워크 variation을 모두 생산하고 이를 이용함으로써, 모델의 구조를 식별할 수 있을 것이다. 그리고 각 모델에 대한 GA결과(아래와 같이 주어지는)를 취합하여 어떤 구조가 가장 좋은 성능을 보이는지를 놓고 결정할 수 있을 것이다.

```
       Avg_Score_Gen        Best_score_Gen
  [1,] "0.183272634920635"  "0.0550556998556999"
  [2,] "0.169272998556999"  "0.0550556998556999"
  [3,] "0.169580124098124"  "0.0423436507936508"
  [4,] "0.148200437229437"  "0.0423436507936508"
  [5,] "0.13649340981241"   "0.0423436507936508"
  [6,] "0.139420212121212"  "0.0423436507936508"
  [7,] "0.131471187590188"  "0.0423436507936508"
  [8,] "0.127832728715729"  "0.0423436507936508"
  [9,] "0.134347287157287"  "0.0423436507936508"
 [10,] "0.135835092352092"  "0.0423436507936508"
 [11,] "0.149200080808081"  "0.0423436507936508"
 [12,] "0.147543392496393"  "0.0296452380952381"
 [13,] "0.139693046176046"  "0.0296452380952381"
 [14,] "0.128549098124098"  "0.0296452380952381"
 ...
```

`+`로 연결된 몇개의 리액션이 보이는데, 이들은 `and`로 연결되어 있다고 설명이 된다.
```
       EGF+TNFa=PI3K Raf+!Akt=Mek Erk+TNFa=Hsp27
EGF                0            0              0
TNFa               0            0              0
Jnk                0            0              0
PI3K               0            0              0
Raf                0            0              0
...

$speciesCompressed
[1] "TRAF6" "p38"   "Ras"

$SplitANDs
$SplitANDs$initialReac
[1] "split1" "split2"

$newANDs
$newANDs$finalReac
[1] "or1" "or2"

$newANDs$`!cJun+TNFa=Jnk`
[1] "!cJun=Jnk" "TNFa=Jnk"

$newANDs$`EGF+TNFa=PI3K`
[1] "EGF=PI3K"  "TNFa=PI3K"

$newANDs$`Raf+!Akt=Mek`
[1] "Raf=Mek"  "!Akt=Mek"

$newANDs$`Erk+TNFa=Hsp27`
[1] "Erk=Hsp27"  "TNFa=Hsp27"

> ToyNCNOcutCompGatexp
```

GA실행결과로써 선택되는 것은 어떤 리액션이 선택되는가와, 그때의 스코어이다.

```
> ToyT1opt
$bString
 [1] 1 1 0 1 1 0 1 1 0 0 0 1 1 1 0 0 0 0

$bScore
0.02963615
```

## Python wrapper

http://www.ebi.ac.uk/~cokelaer/cellnopt/wrapper/quickstart.html#converting-r-script-into-python

```python
from cellnopt import wrapper
# b = CNORbool(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
b = wrapper.CNORbool(wrapper.cnodata("ToyModelMMB2.sif"), wrapper.cnodata("ToyDataMMB.csv"))
b.preprocessing() # compression/expansion/cutNONC species
b.gaBinaryT1(popsize=50, maxgens=10)
b.plotFit()
b.cutAndPlotResultsT1()

# in ipython, you can check bScore with:
In [54]: x = list(tuple(b.T1opt.bScore))
In [55]: x
Out[55]: [0.02963614718614719]
In [56]:
```

# CellNetOpt
