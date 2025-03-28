# 气藏开发软件qcs

## 主要功能

- 数据查询功能
- 常用的气藏及井筒相关计算功能

## 技术栈

- 后端：Python
- 前端：Flask
- 数据库：Mysql

## 用到的第三方库

- python-dotenv
  为了保证重要信息（eg.账号、数据库密码登）不泄露，需要把关键数据写入.env文件中，采用python-dotenv库来加载环境变量中的.env文件
- pymysql
  用于Python与Mysql交互

## Todo List

- [ ] 配置SSL
- [ ] 编写计算函数
