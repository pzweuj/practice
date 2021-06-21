# coding=utf-8
# pzw

import yaml

# 基因型前小后大规则写，如"7/8"
# 三联体是否符合孟德尔定律
def isMendel_triple(alleged_father, mother, child):
    father_type = alleged_father.split("/")
    mother_type = mother.split("/")

    typeList = []
    for f in father_type:
        for m in mother_type:
            typeList.append(f + "/" + m)
            typeList.append(m + "/" + f)

    if child in typeList:
        return True
    else:
        return False

# 二联体是否复核孟德尔定律
def isMendel_double(alleged, child):
    alleged_type = alleged.split("/")
    child_type = child.split("/")
    n = 0
    for c in child_type:
        if c in alleged_type:
            n += 1
    if n != 0:
        return True
    else:
        return False

# 人群频率读取
# freq = frequencyDict("frequency.yml")["D6S1017"]["7"]
def frequencyDict(yamlFile):
    with open(yamlFile, "r") as stream:
        try:
            fqDict = yaml.safe_load(stream)
            return fqDict
        except yaml.YAMLError as exc:
            print(exc)

# 三联体PI计算
def PI_triple(loci, alleged_father, mother, child, frqDB, alleged_gender="male"):
    freqDatabase = frequencyDict(frqDB)
    lociFreq = freqDatabase[loci]

    child_type = child.split("/")
    mother_type = mother.split("/")
    alleged_father_type = alleged_father.split("/")        
    c1 = child_type[0]
    c2 = child_type[1]
    m1 = mother_type[0]
    m2 = mother_type[1]
    af1 = alleged_father_type[0]
    af2 = alleged_father_type[1]

    # 符合孟德尔
    if isMendel_triple(alleged_father, mother, child):
        # 母亲和孩子基因型相同并且等位基因不同，不能确定哪个等位基因来源于父亲
        if mother == child and c1 != c2:
            predict_father = ",".join([c1, c2])
            alleged_father_type.append(c1)
            alleged_father_type.append(c2)
            # 如果存在其他等位基因
            if len(set(alleged_father_type)) > 2:
                pi = 1 / (2 * (lociFreq[c1] + lociFreq[c2]))
            # 不存在其他
            else:
                pi = 2 / (2 * (lociFreq[c1] + lociFreq[c2]))

        # 能确定哪个等位基因来源于父亲
        else:
            if child == mother:
                predict_father = c1
            else:
                predict_father = list(set(child_type).difference(set(mother_type)))[0]            
            # 统计疑似父亲中有几个相同等位基因
            countPF = alleged_father_type.count(predict_father)
            pi = countPF / (2 * lociFreq[predict_father])

    # 不符合孟德尔
    else:
        predict_father = "Nan"
        u_male = 0.0002
        u_female = 0.0005
        if alleged_gender == "female":
            u = u_female
        else:
            u = u_male

        # 判断孩子里是否有来自母亲的
        if c1 in mother_type or c2 in mother_type:
            u = 0.002
            # 判断孩子和母亲基因型是否一致
            # 即两条都有可能来源于父亲
            if child == mother:
                step1 = abs(float(c1) - float(af1))
                step2 = abs(float(c1) - float(af2))
                step3 = abs(float(c2) - float(af1))
                step4 = abs(float(c2) - float(af2))
                step = [step1, step2, step3, step4]
                minStep = min(step)
                minStepCount = step.count(minStep)
                pi = minStepCount * u / (4 * ((lociFreq[c1] + lociFreq[c2]) * (10 ** (minStep - 1))))
            # 不一致表明有且只有一个和母亲相同
            else:
                # 取差集，即other来源于父亲的突变
                other = list(set(child_type).difference(set(mother_type)))[0]
                step1 = abs(float(other) - float(af1))
                step2 = abs(float(other) - float(af2))
                step = [step1, step2]
                minStep = min(step)
                minStepCount = step.count(minStep)
                pi = minStepCount * u / (4 * (lociFreq[other] * (10 ** (minStep - 1))))
        else:
            pi = 0
            print("ERROR，请确认基因型！")

    return [pi, predict_father]


# 二联体PI计算
def PI_double(loci, alleged, child, frqDB, alleged_gender="male"):
    freqDatabase = frequencyDict(frqDB)
    lociFreq = freqDatabase[loci]

    child_type = child.split("/")
    alleged_type = alleged.split("/")  
    c1 = child_type[0]
    c2 = child_type[1]
    at1 = alleged_type[0]
    at2 = alleged_type[1]

    # 符合孟德尔定律
    if isMendel_double(alleged, child):
        # 孩子基因型和疑似相同
        if child == alleged:
            if c1 == c2:
                pi = 2 * (lociFreq[c1] + lociFreq[c2]) / (4 * lociFreq[c1] * lociFreq[c2])
            else:
                pi = (lociFreq[c1] + lociFreq[c2]) / (4 * lociFreq[c1] * lociFreq[c2])
        # 孩子基因型和疑似不同
        else:
            if c1 == c2:
                pi = 1 / (2 * lociFreq[c1])
            else:
                if at1 == at2:
                    pi = 1 / (2 * lociFreq[at1])
                else:
                    other = list(set(child_type).intersection(set(alleged_type)))[0]
                    pi = 1 / (4 * lociFreq[other])
    # 不符合孟德尔
    else:
        u_male = 0.0002
        u_female = 0.0005
        if alleged_gender == "female":
            u = u_female
        else:
            u = u_male

        step1 = abs(float(c1) - float(at1))
        step2 = abs(float(c1) - float(at2))
        step3 = abs(float(c2) - float(at1))
        step4 = abs(float(c2) - float(at2))
        step = [step1, step2, step3, step4]
        minStep = min(step)
        minStepDict = {}
        minStepDict["CHILD1>ALLEGED1"] = step1
        minStepDict["CHILD1>ALLEGED2"] = step2
        minStepDict["CHILD2>ALLEGED1"] = step3
        minStepDict["CHILD2>ALLEGED2"] = step4
        minStepC1 = 0
        minStepC2 = 0
        for i in minStepDict.keys():
            if minStepDict[i] == minStep:
                if "CHILD1" in i:
                    minStepC1 += 1
                if "CHILD2" in i:
                    minStepC2 += 1
        pi_1 = minStepC1 * u / (8 * lociFreq[c1] * (10 ** (minStep - 1)))
        pi_2 = minStepC2 * u / (8 * lociFreq[c2] * (10 ** (minStep - 1)))
        pi = pi_1 + pi_2

    return pi

# 导入结果表格
# #loci  alleged_father    mother  child
# D6S1017   7/8 7/7 8/8
# D1S1627   10/11   12/13   11/12
def calCPI(resultsFile, frqDB, alleged_gender="male"):
    results = open(resultsFile, "r")
    CPI = 1
    for line in resultsFile:
        if not line.startswith("#"):
            lines = line.replace("\n", "").split("\t")
            loci = lines[0]
            alleged_father = lines[1]
            mother = lines[2]
            child = lines[3]

            if mother == "":
                PI = PI_double(loci, alleged_father, child, frqDB, alleged_gender)
            else:
                PI = PI_triple(loci, alleged_father, mother, child, frqDB, alleged_gender)[0]            
            CPI = CPI * PI
    results.close()

    return CPI