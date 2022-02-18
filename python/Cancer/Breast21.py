# coding=utf-8
# pzw
# 20220218
# 乳腺癌21基因RS计算公式

# https://github.com/bhklab/genefu
def scale_genefu(geneDict):
    scaleDict = {}
    geneList = [
        "GRB7", "HER2",
        "ER", "PgR", "BCL2", "SCUBE2",
        "Survivin", "Ki67", "MYBL2", "CCNB1", "STK15",
        "CTSL2", "MMP11",
        "CD68", "GSTM1", "BAG1"
    ]
    geneCT = []
    for g in geneList:
        geneCT.append(geneDict[g])
    for g in geneList:
        scaleCT = (geneDict[g] - min(geneCT)) / (max(geneCT) - min(geneCT)) * 15
        scaleDict[g] = scaleCT
    return scaleDict

# 迪安、美中基因
def scale_d(geneDict):
    referenceCT = [geneDict["ACTB"], geneDict["GAPDH"], geneDict["RPLP0"], geneDict["GUS"], geneDict["TFRC"]]
    meanReferenceCT = sum(referenceCT) / len(referenceCT)
    scaleDict = {}
    for g in geneDict.keys():
        deltaCT = geneDict[g] - meanReferenceCT
        scaleCT = 10 - deltaCT
        scaleDict[g] = scaleCT
    return scaleDict

# test scale
def scale(geneDict):
    referenceCT = [geneDict["ACTB"], geneDict["GAPDH"], geneDict["RPLP0"], geneDict["GUS"], geneDict["TFRC"]]
    meanReferenceCT = sum(referenceCT) / len(referenceCT)
    scaleDict = {}
    for g in geneDict.keys():
        deltaCT = geneDict[g] - meanReferenceCT
        scaleCT = 2 ** (-1 * deltaCT)
        scaleDict[g] = scaleCT
    return scaleDict

def BreastRS(geneDict):
    # GRB7 group score = 0.9 × GRB7 + 0.1 × HER2 (if result < 8, then score is 8)
    HER2_gs = 0.9 * geneDict["GRB7"] + 0.1 * geneDict["HER2"]
    if HER2_gs < 8:
        HER2_gs = 8
    
    # ER group score = (0.8 × ER + 1.2 × PGR + BCL2 + SCUBE2) ÷ 4
    ER_gs = (0.8 * geneDict["ER"] + 1.2 * geneDict["PgR"] + geneDict["BCL2"] + geneDict["SCUBE2"]) / 4
    
    # proliferation group score = (Survivin + KI67 + MYBL2 + CCNB1 + STK15) ÷ 5 (if result < 6.5, then score is 6.5)
    Prolifieration_gs = (geneDict["Survivin"] + geneDict["Ki67"] + geneDict["MYBL2"] + geneDict["CCNB1"] + geneDict["STK15"]) / 5
    if Prolifieration_gs < 6.5:
        Prolifieration_gs = 6.5
    
    # invasion group score = (CTSL2 + MMP11) ÷ 2
    Invasion_gs = (geneDict["CTSL2"] + geneDict["MMP11"]) / 2

    # RSU = + 0.47 × GRB7 group score – 0.34 × ER group score + 1.04 × proliferation group score + 0.10 × invasion group score + 0.05 × CD68 – 0.08 × GSTM1 – 0.07 × BAG1
    RSu = 0.47 * HER2_gs - 0.34 * ER_gs + 1.04 * Prolifieration_gs + 0.1 * Invasion_gs + 0.05 * geneDict["CD68"] - 0.08 * geneDict["GSTM1"] - 0.07 * geneDict["BAG1"]
    
    # RS=0 if RSU<0; RS=20×(RSU–6.7) if 0≤RSU≤100; and RS=100 if RSU>100
    if RSu <= 0:
        RS = 0
    elif RSu >= 100:
        RS = 100
    else:
        RS = 20 * (RSu - 6.7)
    return [RSu, RS]

# Test 9.70 60.00
ctDict = {
    "GRB7": 36.59, "HER2": 33.30,
    "ER": 34.78, "PgR": 34.86, "BCL2": 35.90, "SCUBE2": 37.84,
    "Survivin": 33.82, "Ki67": 35.00, "MYBL2": 34.76, "CCNB1": 33.21, "STK15": 35.42,
    "CTSL2": 38.55, "MMP11": 31.80,
    "CD68": 28.11, "GSTM1": 36.94, "BAG1": 34.83,
    "ACTB": 24.74, "GAPDH": 27.23, "RPLP0": 28.71, "GUS": 36.33, "TFRC": 31.21
}

# Test2 8.42 34.47
# ctDict2 = {
#     "GRB7": 29.27, "HER2": 27.95,
#     "ER": 23.16, "PgR": 26.75, "BCL2": 24.81, "SCUBE2": 30.80,
#     "Survivin": 28.63, "Ki67": 30.90, "MYBL2": 32.92, "CCNB1": 29.59, "STK15": 28.94,
#     "CTSL2": 32.32, "MMP11": 25.96,
#     "CD68": 26.51, "GSTM1": 26.25, "BAG1": 26.68,
#     "ACTB": 20.95, "GAPDH": 22.96, "RPLP0": 24.18, "GUS": 27.08, "TFRC": 28.20
# }

# Test3 7.90 24.09
# ctDict3 = {
#     "GRB7": 26.62, "HER2": 24.42,
#     "ER": 23.54, "PgR": 24.59, "BCL2": 25.46, "SCUBE2": 23.84,
#     "Survivin": 23.84, "Ki67": 26.83, "MYBL2": 30.45, "CCNB1": 26.22, "STK15": 27.74,
#     "CTSL2": 30.56, "MMP11": 19.53,
#     "CD68": 23.52, "GSTM1": 23.68, "BAG1": 26.83,
#     "ACTB": 18.25, "GAPDH": 20.87, "RPLP0": 20.32, "GUS": 22.62, "TFRC": 24.91
# }

s = scale_d(ctDict)
print(s)
rs = BreastRS(s)
print(rs)
