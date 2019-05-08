# pzw
# 20190508

import numpy as np
import pandas as pd

# 导入
dataset = pd.read_csv("Data.csv")
# print(dataset)
X = dataset.iloc[ : , :-1].values
Y = dataset.iloc[ : , 3].values

# 处理丢失数据
from sklearn.impute import SimpleImputer
imputer = SimpleImputer(missing_values=np.NAN, strategy="mean", fill_value=None, verbose=0, copy=True)
imputer = imputer.fit(X[ : , 1:3])
X[ : , 1:3] = imputer.transform(X[ : , 1:3])

# 解析分类数据
from sklearn.preprocessing import LabelEncoder
labelencoder_X = LabelEncoder()
X[ : , 0] = labelencoder_X.fit_transform(X[ : , 0])

# 创建虚拟变量
from sklearn.preprocessing import OneHotEncoder
onehotencoder = OneHotEncoder(categorical_features = [0])
X = onehotencoder.fit_transform(X).toarray()
labelencoder_Y = LabelEncoder()
Y =  labelencoder_Y.fit_transform(Y)

# print(X)

# # 拆分数据集为训练集合和测试集合
from sklearn.model_selection import train_test_split
X_train, X_test, Y_train, Y_test = train_test_split( X , Y , test_size = 0.2, random_state = 0)

# # 特征量化
from sklearn.preprocessing import StandardScaler
sc_X = StandardScaler()
X_train = sc_X.fit_transform(X_train)
X_test = sc_X.transform(X_test)

print(X)