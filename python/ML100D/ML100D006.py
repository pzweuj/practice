# coding=utf-8
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


## 预处理
# 导入
dataset = pd.read_csv("Social_Network_Ads.csv")
# 年龄-收入
X = dataset.iloc[:, [2, 3]].values
# 购买情况
Y = dataset.iloc[:, 4].values

# 将数据集分成训练集和测试集
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.25, random_state=0)

# 特征缩放
from sklearn.preprocessing import StandardScaler
sc = StandardScaler()

X_train = X_train.astype(np.float64)
X_test = X_test.astype(np.float64)

X_train = sc.fit_transform(X_train)
X_test = sc.transform(X_test)

## 逻辑回归模型
# 将逻辑回归应用于训练集
from sklearn.linear_model import LogisticRegression
classifier = LogisticRegression(solver="liblinear")
classifier.fit(X_train, y_train)

## 预测
# 预测测试集结果
y_pred = classifier.predict(X_test)

## 评估预测
# 生成混淆矩阵
from sklearn.metrics import confusion_matrix
cm = confusion_matrix(y_test, y_pred)

# 可视化
from matplotlib.colors import ListedColormap
X_set, y_set = X_train, y_train
X1, X2 = np.meshgrid(np.arange(start=X_set[:,0].min()-1, stop=X_set[:, 0].max()+1, step=0.01),
	np.arange(start=X_set[:, 1].min()-1, stop=X_set[:, 1].max()+1, step=0.01))
plt.contourf(X1, X2, classifier.predict(np.array([X1.ravel(), X2.ravel()]).T).reshape(X1.shape),
	alpha = 0.75, cmap = ListedColormap(("red", "green")))
plt.xlim(X1.min(), X1.max())
plt.ylim(X2.min(), X2.max())
for i, j in enumerate(np.unique(y_set)):
	plt.scatter(X_set[y_set==j,0], X_set[y_set==j,1],
		c = ListedColormap(("red", "green"))(i), label=j)

plt.title(" LOGISTIC(Training set)")
plt.xlabel(" Age")
plt.ylabel(" Estimated Salary")
plt.legend()
plt.show()

X_set, y_set = X_test, y_test
X1, X2=np.meshgrid(np.arange(start=X_set[:, 0].min()-1, stop=X_set[:, 0].max()+1, step=0.01),
	np.arange(start=X_set[:, 1].min()-1, stop=X_set[:, 1].max()+1, step=0.01))

plt.contourf(X1, X2, classifier.predict(np.array([X1.ravel(), X2.ravel()]).T).reshape(X1.shape),
	alpha=0.75, cmap=ListedColormap(("red", "green")))
plt.xlim(X1.min(), X1.max())
plt.ylim(X2.min(), X2.max())
for i, j in enumerate(np.unique(y_set)):
	plt.scatter(X_set[y_set==j, 0], X_set[y_set==j, 1],
		c=ListedColormap(("red", "green"))(i), label=j)

plt.title(" LOGISTIC(Test set)")
plt.xlabel(" Age")
plt.ylabel(" Estimated Salary")
plt.legend()
plt.show()