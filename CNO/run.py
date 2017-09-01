from cellnopt import wrapper
# b = CNORbool(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
b = wrapper.CNORbool(wrapper.cnodata("ToyModelMMB2.sif"),
        wrapper.cnodata("ToyDataMMB.csv"))
b.preprocessing() # compression/expansion/cutNONC species
b.gaBinaryT1(popsize=50, maxgens=10)


