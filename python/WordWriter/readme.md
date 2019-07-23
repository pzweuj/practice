#### 依赖
```bash
pip install python-docx
pip install pandas
```


#### 使用方法
```python
## python2
import WordWriter

## python3
import WordWriter3
```

在模板word中创建标签，形式为**#[xxxx]#**。

其中表格标签必须是#[TABLE-xxx]#，

表格内容标签必须是#[TBS-xxx]#，

图片标签必须是#[IMAGE-xxx]#，支持定义图片大小#[IMAGE-xxx-(30,40)]#，

页眉标签必须是#[HEADER-xxx]#，页脚标签必须是#[FOOTER-xxx]#，

其他文本内容标签可自定义#[xxxx]#。

#### 实例
test.docx是自定义模板。

```python
# python2
import WordWriter

# 定义一个字典放入所有内容
testDict = {}
testDict["#[HEADER-1]#"] = "模板测试"
testDict["#[HEADER-2]#"] = "2019年7月18日"
testDict["#[NAME]#"] = "测试模板"
testDict["#[fullParagraph]#"] = "这是一段测试段落，通过WordWriter输入。"
testDict["#[TBS-1]#"] = "未突变"
testDict["#[FOOTER]#"] = "页脚测试"

# 表格标签，图片标签指定的是文件，其中表格标签输入文件是以tab分割的文本，无标题
testDict["#[TABLE-1]#"] = "testTable.txt"
testDict["#[IMAGE-1-(30,30)]#"] = "testPicture.png"
testDict["#[IMAGE-2]#"] = "testPicture.png"

# 使用主函数进行报告填充
WordWriter.WordWriter("test.docx", "testOut.docx", testDict)
```

python3

```python
# python3
import WordWriter3

# 定义一个字典放入所有内容
testDict = {}
testDict["#[HEADER-1]#"] = "模板测试"
testDict["#[HEADER-2]#"] = "2019年7月18日"
testDict["#[NAME]#"] = "测试模板"
testDict["#[fullParagraph]#"] = "这是一段测试段落，通过WordWriter输入。"
testDict["#[TBS-1]#"] = "未突变"
testDict["#[FOOTER]#"] = "页脚测试"

# 表格标签，图片标签指定的是文件，其中表格标签输入文件是以tab分割的文本，无标题
testDict["#[TABLE-1]#"] = "testTable.txt"
testDict["#[IMAGE-1-(30,30)]#"] = "testPicture.png"
testDict["#[IMAGE-2]#"] = "testPicture.png"

# 使用主函数进行报告填充
WordWriter3.WordWriter("test.docx", "testOut.docx", testDict)
```