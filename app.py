import os
from dotenv import load_dotenv
from flask import Flask, render_template, request, redirect, url_for,jsonify,json
import pymysql
from pymysql.cursors import DictCursor
import fuction
from flask_login import LoginManager, UserMixin, login_user, login_required

# 加载环境变量
load_dotenv()

# 创建Flask实例
app = Flask(__name__)
app.secret_key = os.getenv('FLASK_SECRET_KEY')

# 配置Flask-Login
login_manager = LoginManager()
login_manager.init_app(app)
login_manager.login_view = 'login_page'  # 设置登录页面的endpoint

# 用户类
class User(UserMixin):
    def __init__(self, id, username):
        self.id = id
        self.username = username

# 用户加载函数
@login_manager.user_loader
def load_user(user_id):
    try:
        connection = pymysql.connect(**get_db_config())
        with connection.cursor() as cursor:
            sql = "SELECT * FROM users WHERE id = %s"
            cursor.execute(sql, (user_id,))
            user = cursor.fetchone()
            if user:
                return User(user['id'], user['username'])
            return None
    except pymysql.MySQLError as e:
        app.logger.error(f"数据库错误: {str(e)}")
        return None
    finally:
        if 'connection' in locals():
            connection.close()

def get_db_config():
    """安全获取数据库配置"""
    required_keys = ['DB_HOST', 'DB_USER', 'DB_PASSWORD', 'DB_NAME']
    missing = [key for key in required_keys if not os.getenv(key)]
    if missing:
        raise EnvironmentError(f"缺少必需的环境变量: {', '.join(missing)}")

    config = {
        'host': os.getenv('DB_HOST'),
        'user': os.getenv('DB_USER'),
        'password': os.getenv('DB_PASSWORD'),
        'db': os.getenv('DB_NAME'),
        'charset': 'utf8mb4',
        'cursorclass': pymysql.cursors.DictCursor
    }
    return config

@app.route('/')
def login_page():
    return render_template('login.html')

@app.route('/login', methods=['POST'])
def login():
    username = request.form['username']
    password = request.form['password']

    try:
        connection = pymysql.connect(**get_db_config())
        with connection.cursor() as cursor:
            sql = "SELECT * FROM users WHERE username = %s AND password = %s"
            cursor.execute(sql, (username, password))
            user = cursor.fetchone()

        if user:
            # 创建用户对象并登录
            user_obj = User(user['id'], user['username'])
            login_user(user_obj)
            return redirect(url_for('main_cal'))
        else:
            return render_template('login.html', error='用户名或密码错误')

    except pymysql.MySQLError as e:
        app.logger.error(f"数据库错误: {str(e)}")
        return render_template('login.html', error='数据库错误')
    except EnvironmentError as e:
        app.logger.error(f"配置错误: {str(e)}")
        return render_template('login.html', error='服务器配置错误')
    finally:
        if 'connection' in locals():
            connection.close()

@app.route('/main_cal')
@login_required  # 使用Flask-Login的登录验证装饰器
def main_cal():
    return render_template('main_cal.html')

@app.route('/success', methods=['POST', 'GET'])
@login_required  # 使用Flask-Login的登录验证装饰器
def success():
    if request.method == 'POST':
        well_name = request.form.get('well_input')
        try:
            connection = pymysql.connect(**get_db_config())
            with connection.cursor() as cursor:
                sql = "SELECT * FROM dadiaocha WHERE wellname = %s"
                cursor.execute(sql, (well_name,))
                result = cursor.fetchall()
                # ##把结果重新构建一个列表，以显示更直观
                # result = result[0]
                # result = list(result.values())
                # # 假设fuction.z是一个计算函数
                # result = fuction.z(result[18], result[19], result[16], 5)

            # 确保返回列表格式（即使空结果）
            return render_template('success.html',results=result,searched_well=well_name,error=None)

        except pymysql.MySQLError as e:
            return render_template('success.html',results=[],error=str(e))
        finally:
            if 'connection' in locals():
                connection.close()

    elif request.method == 'GET':
        return render_template('success.html',results=[],error=None)

@app.route('/lookfor_month_data', methods=['GET'])
def lookfor_month_data_page():
    """返回包含Handsontable的前端页面"""
    return render_template('lookfor_month_data.html')

@app.route('/api/search', methods=['POST'])
def search_well():
    """API端点，返回JSON格式的查询结果（包含所有可能字段）"""
    well_name = request.form.get('wellname')
    if not well_name:
        return jsonify({
            "success": False,
            "data": [],
            "error": "Well name is required"
        }), 400
    
    try:
        connection = pymysql.connect(**get_db_config())
        with connection.cursor() as cursor:
            # 1. 先查询符合条件的数据
            sql = "SELECT * FROM monthlypro WHERE JH = %s"
            cursor.execute(sql, (well_name,))
            result = cursor.fetchall()
            
            # 2. 获取表的所有字段（确保即使某些记录缺少字段也能显示表头）
            cursor.execute("DESCRIBE monthlypro")
            all_fields = [field['Field'] for field in cursor.fetchall()]
            
            # 3. 规范化数据（确保每条记录包含所有字段）
            normalized_data = []
            for row in result:
                normalized_row = {field: row.get(field, None) for field in all_fields}
                normalized_data.append(normalized_row)
            
            return jsonify({
                "success": True,
                "data": normalized_data,
                "columns": all_fields,  # 返回所有字段名用于前端表头
                "error": None
            })
            
    except pymysql.MySQLError as e:
        return jsonify({
            "success": False,
            "data": [],
            "error": str(e)
        }), 500
    finally:
        if 'connection' in locals():
            connection.close()


if __name__ == '__main__':
    app.run(debug=True,host='0.0.0.0')