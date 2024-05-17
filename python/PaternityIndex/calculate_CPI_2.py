# coding=utf-8
# pzw
# 20240517
# GPT优化代码，未验证

import yaml

# 基因型前小后大规则写，如"7/8"
# 三联体是否符合孟德尔定律
def isMendel_triple(alleged_father, mother, child):
    father_alleles = set(alleged_father.split("/"))
    mother_alleles = set(mother.split("/"))
    
    # 生成所有可能的组合
    possible_combinations = {f + "/" + m for f in father_alleles for m in mother_alleles} | \
                            {m + "/" + f for f in father_alleles for m in mother_alleles}
    
    return child in possible_combinations

# 二联体是否复核孟德尔定律
def isMendel_double(alleged, child):
    alleged_alleles = set(alleged.split("/"))
    child_alleles = set(child.split("/"))
    
    # 使用集合的交集操作来查找共享的等位基因
    shared_alleles = alleged_alleles.intersection(child_alleles)
    
    # 使用 any 函数检查是否存在共享的等位基因
    return any(allele in shared_alleles for allele in child_alleles)

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
    c1, c2 = child_type
    m1, m2 = mother_type
    af1, af2 = alleged_father_type

    # 符合孟德尔
    if isMendel_triple(alleged_father, mother, child):
        # 母亲和孩子基因型相同并且等位基因不同，不能确定哪个等位基因来源于父亲
        if mother == child and c1 != c2:
            predict_father = ",".join([c1, c2])
            alleged_father_type.extend([c1, c2])
            pi = 1 / (2 * (lociFreq[c1] + lociFreq[c2])) if len(set(alleged_father_type)) > 2 else 2 / (2 * (lociFreq[c1] + lociFreq[c2]))

        # 能确定哪个等位基因来源于父亲
        else:
            predict_father = c1 if child == mother else next(iter(set(child_type).difference(set(mother_type))))
            countPF = alleged_father_type.count(predict_father)
            pi = countPF / (2 * lociFreq[predict_father])

    # 不符合孟德尔
    else:
        predict_father = "Nan"
        u_male, u_female = 0.0002, 0.0005
        u = u_female if alleged_gender == "female" else u_male

        # 判断孩子里是否有来自母亲的
        if c1 in mother_type or c2 in mother_type:
            u, other = 0.002, next(iter(set(child_type).difference(set(mother_type))))

            # 判断孩子和母亲基因型是否一致
            if child == mother:
                step = [abs(float(c1) - float(af1)), abs(float(c1) - float(af2)), abs(float(c2) - float(af1)), abs(float(c2) - float(af2))]
                minStep = min(step)
                minStepCount = step.count(minStep)
                pi = minStepCount * u / (4 * ((lociFreq[c1] + lociFreq[c2]) * (10 ** (minStep - 1))))
            else:
                step = [abs(float(other) - float(af1)), abs(float(other) - float(af2))]
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
    c1, c2 = child_type
    at1, at2 = alleged_type

    # 符合孟德尔定律
    if isMendel_double(alleged, child):
        if child == alleged:
            pi = 2 * (lociFreq[c1] + lociFreq[c2]) / (4 * lociFreq[c1] * lociFreq[c2]) if c1 == c2 else \
                 (lociFreq[c1] + lociFreq[c2]) / (4 * lociFreq[c1] * lociFreq[c2])
        else:
            if c1 == c2:
                pi = 1 / (2 * lociFreq[c1])
            else:
                pi = 1 / (2 * lociFreq[at1]) if at1 == at2 else 1 / (4 * lociFreq[list(set(child_type).intersection(set(alleged_type)))[0]])
    # 不符合孟德尔
    else:
        u_male, u_female = 0.0002, 0.0005
        u = u_female if alleged_gender == "female" else u_male

        step1 = abs(float(c1) - float(at1))
        step2 = abs(float(c1) - float(at2))
        step3 = abs(float(c2) - float(at1))
        step4 = abs(float(c2) - float(at2))
        minStep = min(step1, step2, step3, step4)
        minStepCount = sum(step == minStep for step in [step1, step2, step3, step4])

        pi = minStepCount * u / (8 * lociFreq[c1 if step1 == minStep or step2 == minStep else c2] * (10 ** (minStep - 1)))

    return pi

# 导入结果表格
# #loci  alleged_father    mother  child
# D6S1017   7/8 7/7 8/8
# D1S1627   10/11   12/13   11/12
def calCPI(resultsFile, frqDB, alleged_gender="male"):
    CPI = 1
    with open(resultsFile, "r") as results:
        for line in results:
            if not line.startswith("#"):
                loci, alleged_father, mother, child = line.strip().split("\t")

                if mother == "":
                    PI = PI_double(loci, alleged_father, child, frqDB, alleged_gender)
                else:
                    PI = PI_triple(loci, alleged_father, mother, child, frqDB, alleged_gender)[0]
                    
                CPI *= PI

    return CPI

# end
