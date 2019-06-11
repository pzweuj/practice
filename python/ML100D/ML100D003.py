# coding=utf-8
import pandas as pd
import numpy as np

## 预处理
# 导入
dataset = pd.read_csv('50_Startups.csv')
X = dataset.iloc[ : , :-1].values
Y = dataset.iloc[ : ,  4 ].values
# print X, Y

# 将类别虚拟化
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder



from sklearn.preprocessing import LabelEncoder, OneHotEncoder
labelencoder = LabelEncoder()
X[: , 3] = labelencoder.fit_transform(X[ : , 3])
onehotencoder = OneHotEncoder(categorical_features = [3])
X = onehotencoder.fit_transform(X).toarray()


# 
# from sklearn.preprocessing import LabelEncoder, OneHotEncoder
# labelencoder = LabelEncoder()
# X[: , 3] = labelencoder.fit_transform(X[ : , 3])
# onehotencoder = OneHotEncoder(categories="auto")
# X = onehotencoder.fit_transform(X).toarray()
# print X

# # 躲避陷阱
X = X[: , 1:]
# print X

# 拆分数据集为训练集和测试集
from sklearn.model_selection import train_test_split
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size = 0.2, random_state = 0)
# print X_train, Y_train

## 在训练集上训练多元线性回归模型
from sklearn.linear_model import LinearRegression
regressor = LinearRegression()
regressor.fit(X_train, Y_train)

## 在测试集上预测结果
y_pred = regressor.predict(X_test)
print y_pred