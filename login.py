#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
登录网站并获取登录后的页面内容

功能:
- 使用用户名和密码登录网站
- 获取登录后的HTML页面内容
- 保存登录后的HTML到文件

运行方法:
1. 如果使用 conda 环境: python login.py 或 python3 login.py
2. 如果使用系统 Python: 需要先安装 requests: pip install requests
"""

import time
import re
from urllib.parse import urljoin, urlparse

# 尝试导入 requests，如果失败则提示
try:
    import requests
    REQUESTS_AVAILABLE = True
except ImportError:
    REQUESTS_AVAILABLE = False
    print("警告: requests 模块未安装。网络测试功能将不可用。")
    print("\n安装方法:")
    print("  1. 如果使用 conda 环境: conda activate base 然后 pip install requests")
    print("  2. 如果使用系统 Python: python3 -m pip install requests")
    print("  3. 或者直接使用 conda 环境中的 Python 运行此脚本\n")

def login_and_get_content(login_url, username, password, timeout=10, save_to_file=True, api_url=None):
    """
    登录网站并获取登录后的页面内容
    
    参数:
        login_url: 登录页面URL
        username: 用户名
        password: 密码
        timeout: 超时时间（秒）
        save_to_file: 是否保存HTML到文件（默认True）
        api_url: 可选的登录API端点URL（如果已知，可以手动指定）
    
    返回:
        dict: 包含登录状态和页面内容的字典
    """
    if not REQUESTS_AVAILABLE:
        result = {
            'success': False,
            'status_code': None,
            'content': None,
            'error': 'requests 模块未安装'
        }
        print("\n✗ requests 模块未安装，无法进行登录")
        print("\n安装方法:")
        print("  1. 如果使用 conda 环境: conda activate base 然后 pip install requests")
        print("  2. 如果使用系统 Python: python3 -m pip install requests")
        return result
    
    # 创建session以保持cookies
    session = requests.Session()
    
    # 设置请求头，模拟浏览器访问
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
        'Accept-Language': 'zh-CN,zh;q=0.9,en;q=0.8',
        'Accept-Encoding': 'gzip, deflate, br',
        'Connection': 'keep-alive',
        'Upgrade-Insecure-Requests': '1'
    }
    session.headers.update(headers)
    
    result = {
        'success': False,
        'status_code': None,
        'content': None,
        'encoding': None,
        'final_url': None,
        'error': None
    }
    
    try:
        print(f"\n正在访问登录页面: {login_url}")
        print("-" * 50)
        
        # 第一步：先访问首页，获取必要的cookies（如GDSSID等）
        print(f"  步骤0: 访问首页以获取必要的cookies...")
        base_url = f"{urlparse(login_url).scheme}://{urlparse(login_url).netloc}"
        try:
            home_response = session.get(base_url, timeout=timeout)
            print(f"  ✓ 首页访问成功，状态码: {home_response.status_code}")
        except:
            print(f"  ⚠ 首页访问失败，继续访问登录页面...")
        
        # 第二步：访问登录页面，获取登录表单
        start_time = time.time()
        login_page_response = session.get(login_url, timeout=timeout)
        login_page_time = round((time.time() - start_time) * 1000, 2)
        
        print(f"✓ 登录页面加载成功")
        print(f"  状态码: {login_page_response.status_code}")
        print(f"  响应时间: {login_page_time} ms")
        
        # 检查并显示cookies（可能包含GDSSID、PHPSESSIONID等）
        cookies_dict = session.cookies.get_dict()
        if cookies_dict:
            print(f"  获取到的Cookies: {list(cookies_dict.keys())}")
            # 显示重要的cookies
            important_cookies = ['GDSSID', 'PHPSESSIONID', 'SESSIONID', 'sessionid', 'acw_tc']
            for cookie_name in important_cookies:
                if cookie_name in cookies_dict:
                    cookie_value = cookies_dict[cookie_name]
                    print(f"    {cookie_name}: {cookie_value[:50]}...")
        
        # 解析登录页面
        login_page_content = login_page_response.text
        
        # 由于这是JavaScript渲染的单页应用，需要从页面中提取真正的登录API端点
        # 如果用户手动指定了API端点，直接使用
        login_api_url = api_url
        if login_api_url:
            print(f"  使用手动指定的API端点: {login_api_url}")
        else:
            # 方法1: 查找JavaScript代码中的API调用
            login_api_url = None
        
        # 查找常见的API调用模式
        # 查找 fetch, axios, XMLHttpRequest 等网络请求
        api_patterns = [
            r'["\']([^"\']*api[^"\']*login[^"\']*)["\']',  # 包含api和login的字符串
            r'["\']([^"\']*user[^"\']*login[^"\']*)["\']',  # 包含user和login的字符串
            r'["\']([^"\']*auth[^"\']*login[^"\']*)["\']',  # 包含auth和login的字符串
            r'url\s*[:=]\s*["\']([^"\']*login[^"\']*)["\']',  # url: "...login..."
            r'endpoint\s*[:=]\s*["\']([^"\']*login[^"\']*)["\']',  # endpoint: "...login..."
        ]
        
        found_endpoints = set()
        for pattern in api_patterns:
            matches = re.findall(pattern, login_page_content, re.IGNORECASE)
            for match in matches:
                # 清理匹配结果
                match = match.strip()
                if match and ('login' in match.lower() or 'auth' in match.lower()):
                    # 如果是相对路径，补全
                    if match.startswith('/'):
                        found_endpoints.add(urljoin(base_url, match))
                    elif match.startswith('http'):
                        found_endpoints.add(match)
                    else:
                        found_endpoints.add(urljoin(base_url, '/' + match))
        
        if found_endpoints:
            login_api_url = list(found_endpoints)[0]
            print(f"  在页面JavaScript中找到API端点: {login_api_url}")
            if len(found_endpoints) > 1:
                print(f"  找到多个可能的端点: {list(found_endpoints)}")
        
        # 方法2: 如果没有找到，尝试常见的API路径
        if not login_api_url:
            possible_api_paths = [
                '/api/user/login',
                '/api/user/account/login',
                '/api/auth/login',
                '/api/login',
                '/api/account/login',
                '/user/login',
                '/auth/login',
                '/v2/api/user/login',
                '/v2/api/auth/login'
            ]
            
            for api_path in possible_api_paths:
                if api_path in login_page_content:
                    login_api_url = urljoin(base_url, api_path)
                    print(f"  在页面内容中找到可能的API端点: {login_api_url}")
                    break
        
        # 方法3: 使用已知的正确登录API端点（账户密码登录）
        # 直接使用已知的API端点，就像用户在页面上选择"账户密码登录"后点击登录一样
        if not login_api_url:
            login_api_url = "https://apigateway.gaodun.com/api/v1/ucenter/login"
            print(f"  使用账户密码登录API端点: {login_api_url}")
        
        # 确定要使用的API端点
        if api_url:
            api_urls_to_try = [api_url]
            print(f"  使用手动指定的API端点: {api_url}")
        else:
            # 直接使用已知的正确端点（模拟用户在页面上选择账户密码登录）
            api_urls_to_try = [login_api_url]
            print(f"  模拟页面操作：选择账户密码登录，提交到: {login_api_url}")
        
        # 准备登录数据（根据实际请求参数）
        # 根据提供的实际参数，登录使用form-data格式，字段名为 'user' 和 'password'
        from urllib.parse import quote
        
        # 基础登录数据（form-data格式）
        login_data_form = {
            'user': username,
            'password': password,
            'appid': '210804',
            'operateSource': 'PC',
            'deviceType': 'PC',
            'refer': login_url,
            'languageType': 'CH'
        }
        
        # 尝试从页面或cookies中获取其他必要的参数
        # GDSSID, PHPSESSIONID, webUmidToken, uaToken 等可能需要从页面或cookies中获取
        
        # 从cookies中获取必要的参数（在访问登录页面后，cookies已经设置）
        # cookies_dict已经在上面获取过了，这里直接使用
        if 'GDSSID' in cookies_dict:
            login_data_form['GDSSID'] = cookies_dict['GDSSID']
            print(f"  从Cookies获取GDSSID: {cookies_dict['GDSSID'][:30]}...")
        if 'PHPSESSIONID' in cookies_dict:
            login_data_form['PHPSESSIONID'] = cookies_dict['PHPSESSIONID']
            print(f"  从Cookies获取PHPSESSIONID: {cookies_dict['PHPSESSIONID'][:30]}...")
        
        # 注意：webUmidToken和uaToken可能需要从页面JavaScript中动态生成
        # 如果登录失败，可能需要从页面中提取这些token
        
        # 尝试从页面中提取token等信息
        # 查找可能的token或session ID（更广泛的搜索模式）
        token_patterns = [
            (r'GDSSID["\']?\s*[:=]\s*["\']([^"\']+)["\']', 'GDSSID'),
            (r'PHPSESSIONID["\']?\s*[:=]\s*["\']([^"\']+)["\']', 'PHPSESSIONID'),
            (r'webUmidToken["\']?\s*[:=]\s*["\']([^"\']+)["\']', 'webUmidToken'),
            (r'uaToken["\']?\s*[:=]\s*["\']([^"\']+)["\']', 'uaToken'),
            # 也尝试查找在JavaScript变量中的值
            (r'GDSSID\s*=\s*["\']([^"\']+)["\']', 'GDSSID'),
            (r'PHPSESSIONID\s*=\s*["\']([^"\']+)["\']', 'PHPSESSIONID'),
        ]
        
        for pattern, token_name in token_patterns:
            matches = re.findall(pattern, login_page_content, re.IGNORECASE)
            if matches and token_name not in login_data_form:
                login_data_form[token_name] = matches[0]
                print(f"  从页面提取{token_name}: {matches[0][:30]}...")
        
        # 如果仍然缺少关键参数，尝试生成或使用默认值
        # 注意：某些token可能是客户端JavaScript生成的，无法从页面直接获取
        if 'GDSSID' not in login_data_form:
            print(f"  ⚠ 警告: 未找到GDSSID，可能影响登录")
        if 'PHPSESSIONID' not in login_data_form:
            print(f"  ⚠ 警告: 未找到PHPSESSIONID，可能影响登录")
        
        # 构建JSON格式的数据（备用）
        login_data_json = {
            'user': username,
            'password': password,
            'appid': '210804',
            'operateSource': 'PC',
            'deviceType': 'PC',
            'refer': login_url,
            'languageType': 'CH'
        }
        
        print(f"\n正在模拟页面操作：账户密码登录...")
        print(f"  步骤1: 已访问登录页面 ✓")
        print(f"  步骤2: 选择账户密码登录方式 ✓")
        print(f"  步骤3: 填写用户名和密码...")
        print(f"    用户名字段: user")
        print(f"    密码字段: password")
        print(f"    用户名: {username}")
        print(f"  步骤4: 准备提交登录表单...")
        print(f"    登录数据字段: {list(login_data_form.keys())}")
        
        # 确定要尝试的API端点列表
        # 已知的正确登录API端点
        correct_api_url = "https://apigateway.gaodun.com/api/v1/ucenter/login"
        
        if api_url:
            # 用户手动指定了API端点
            api_urls_to_try = [api_url]
            print(f"  使用手动指定的API端点: {api_url}")
        elif login_api_url == correct_api_url:
            # 使用已知的正确端点
            api_urls_to_try = [correct_api_url]
            print(f"  使用已知的正确API端点: {correct_api_url}")
        elif login_api_url and login_api_url != urljoin(base_url, '/api/user/login'):
            # 从页面中找到了API端点
            print(f"  使用从页面中找到的API端点: {login_api_url}")
            api_urls_to_try = [login_api_url]
        else:
            # 优先使用已知的正确端点，然后尝试其他常见端点
            print(f"  将尝试多个API端点（优先使用已知的正确端点）")
            api_urls_to_try = [
                correct_api_url,  # 优先使用已知的正确端点
                urljoin(base_url, '/api/user/login'),
                urljoin(base_url, '/api/user/account/login'),
                urljoin(base_url, '/api/auth/login'),
                urljoin(base_url, '/api/login'),
                urljoin(base_url, '/v2/api/user/login')
            ]
        
        # 更新请求头
        session.headers.update({
            'Referer': login_url,
            'Origin': base_url
        })
        
        # 第二步：尝试多个API端点和数据格式
        login_response = None
        login_time = 0
        
        # 直接使用form-data格式（模拟页面表单提交）
        # 就像用户在页面上选择"账户密码登录"后点击登录按钮一样
        data_formats = [
            ('form', login_data_form)  # 使用form-data格式，模拟表单提交
        ]
        
        login_success = False
        for api_url in api_urls_to_try:
            if login_success:
                break
            for format_type, data in data_formats:
                try:
                    print(f"  提交登录请求到: {api_url} (格式: {format_type})")
                    
                    # 显示实际提交的数据（隐藏密码）
                    data_display = {k: (v[:20] + '...' if len(str(v)) > 20 else v) for k, v in data.items()}
                    if 'password' in data_display:
                        data_display['password'] = '***'
                    print(f"  提交的数据字段: {list(data.keys())}")
                    print(f"  数据预览: {data_display}")
                    
                    start_time = time.time()
                    
                    if format_type == 'json':
                        session.headers.update({
                            'Content-Type': 'application/json',
                            'Accept': 'application/json, text/plain, */*',
                        })
                        response = session.post(api_url, json=data, timeout=timeout, allow_redirects=False)
                    else:  # form格式
                        session.headers.update({
                            'Content-Type': 'application/x-www-form-urlencoded',
                            'Accept': 'application/json, text/html, */*',
                        })
                        response = session.post(api_url, data=data, timeout=timeout, allow_redirects=False)
                    
                    login_time = round((time.time() - start_time) * 1000, 2)
                    
                    # 输出响应信息用于调试
                    response_preview = response.text[:300] if hasattr(response, 'text') else str(response.content[:300])
                    print(f"    状态码: {response.status_code}, 响应预览: {response_preview}")
                    
                    # 检查响应是否成功
                    if response.status_code in [200, 201, 302]:
                        try:
                            resp_json = response.json()
                            print(f"    JSON响应: {resp_json}")
                            # 检查是否登录成功
                            if isinstance(resp_json, dict):
                                # 检查各种成功标识
                                status = resp_json.get('status')
                                code = resp_json.get('code')
                                info = resp_json.get('info', '')
                                result_msg = resp_json.get('result', '')
                                
                                # 如果返回"参数错误"，说明API端点正确但参数有问题
                                if '参数错误' in str(info) or '参数错误' in str(result_msg):
                                    print(f"    ⚠ API返回参数错误，可能缺少必要的参数")
                                    print(f"    当前提交的参数: {list(data.keys())}")
                                    print(f"    可能缺少的参数: GDSSID, PHPSESSIONID, webUmidToken, uaToken")
                                    # 保存响应以便后续分析
                                    login_response = response
                                    login_api_url = api_url
                                    # 不设置login_success为True，继续尝试或返回错误
                                
                                # 检查成功标识
                                elif code in [200, 0] or status == 200 or resp_json.get('success') == True or 'accessToken' in resp_json or 'token' in resp_json or 'data' in resp_json:
                                    login_response = response
                                    login_api_url = api_url
                                    login_success = True
                                    print(f"  ✓ 登录成功！使用端点: {api_url}")
                                    break
                        except Exception as json_err:
                            # 如果不是JSON，检查状态码和内容
                            if response.status_code == 200:
                                # 检查响应内容中是否有成功标识
                                response_text_lower = response.text.lower() if hasattr(response, 'text') else ""
                                if 'success' in response_text_lower or 'token' in response_text_lower or response.status_code == 302:
                                    login_response = response
                                    login_api_url = api_url
                                    login_success = True
                                    print(f"  ✓ 登录成功！使用端点: {api_url} (非JSON响应)")
                                    break
                    
                except requests.exceptions.RequestException as e:
                    continue
                except Exception as e:
                    continue
        
        if not login_response:
            # 如果所有尝试都失败
            result['error'] = "无法找到有效的登录API端点或参数错误"
            print(f"\n✗ 所有登录尝试都失败")
            print(f"\n可能的原因:")
            print(f"  1. API端点不正确")
            print(f"  2. 缺少必要的参数（如GDSSID、PHPSESSIONID、webUmidToken、uaToken等）")
            print(f"  3. 这些token可能需要从页面JavaScript中动态生成")
            print(f"\n建议:")
            print(f"  1. 打开浏览器，访问 {login_url}")
            print(f"  2. 按F12打开开发者工具，切换到Network标签")
            print(f"  3. 在登录页面选择'账户密码登录'，输入用户名和密码，点击登录")
            print(f"  4. 在Network标签中找到登录请求（POST请求）")
            print(f"  5. 查看Request Payload/Form Data，确认所有必需的参数")
            print(f"  6. 查看Request Headers，确认cookies和headers")
            return result
        elif not login_success:
            # 如果API返回了响应但登录失败（如参数错误）
            result['error'] = "登录API返回参数错误，可能缺少必要的token参数"
            print(f"\n⚠ 登录失败: {result['error']}")
            print(f"\n这些token（GDSSID、PHPSESSIONID、webUmidToken、uaToken）可能需要:")
            print(f"  1. 从页面JavaScript中动态生成")
            print(f"  2. 通过特定的API调用获取")
            print(f"  3. 或者这些参数在某些情况下是可选的")
            print(f"\n建议使用浏览器开发者工具查看实际的登录请求，确认所有必需的参数")
            # 仍然返回响应内容，以便用户查看
            if login_response:
                try:
                    result['content'] = login_response.text if hasattr(login_response, 'text') else str(login_response.content)
                except:
                    pass
            return result
        
        result['status_code'] = login_response.status_code
        result['response_time'] = login_time
        
        # 检查登录API响应
        print(f"\n登录API响应:")
        print(f"  状态码: {login_response.status_code}")
        print(f"  响应时间: {login_time} ms")
        
        # 尝试解析JSON响应
        login_api_success = False
        try:
            login_json = login_response.json()
            print(f"  API响应: {login_json}")
            
            # 检查常见的成功标识
            if login_response.status_code == 200:
                # 检查响应中是否包含成功标识
                if isinstance(login_json, dict):
                    if login_json.get('code') == 200 or login_json.get('success') == True or 'token' in login_json or 'data' in login_json:
                        login_api_success = True
                    elif login_json.get('code') == 0 or login_json.get('status') == 'success':
                        login_api_success = True
                else:
                    login_api_success = True
        except:
            # 如果不是JSON，检查响应内容
            response_text = login_response.text[:500]
            print(f"  响应内容预览: {response_text}")
            if login_response.status_code == 200 and ('success' in response_text.lower() or 'token' in response_text.lower()):
                login_api_success = True
        
        if not login_api_success:
            result['error'] = f"登录API返回失败，状态码: {login_response.status_code}"
            print(f"✗ 登录失败: {result['error']}")
            result['content'] = login_response.text if hasattr(login_response, 'text') else str(login_response.content)
            return result
        
        print(f"✓ 登录API调用成功！")
        
        # 登录成功后，访问登录后的页面（通常是首页或用户中心）
        # 恢复Content-Type为HTML
        session.headers.update({
            'Content-Type': 'text/html',
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8'
        })
        
        # 尝试访问登录后的页面（根据日志，登录后会跳转到 /space/）
        possible_after_login_urls = [
            urljoin(base_url, '/space/'),  # 学习空间 - 登录后的主页面
            urljoin(base_url, '/space'),    # 不带斜杠的版本
            urljoin(base_url, '/'),
            urljoin(base_url, '/home'),
            urljoin(base_url, '/dashboard'),
            urljoin(base_url, '/user'),
            urljoin(base_url, '/user/center'),
            urljoin(base_url, '/account')
        ]
        
        after_login_url = None
        after_login_content = None
        
        print(f"\n正在获取登录后的页面...")
        for test_url in possible_after_login_urls:
            try:
                print(f"  尝试访问: {test_url}")
                page_response = session.get(test_url, timeout=timeout)
                if page_response.status_code == 200:
                    after_login_url = test_url
                    after_login_content = page_response.text
                    print(f"  ✓ 成功获取页面: {test_url}")
                    break
            except:
                continue
        
        # 如果没有找到，使用重定向后的URL
        if not after_login_content and login_response.url != login_api_url:
            try:
                after_login_url = login_response.url
                after_login_content = login_response.text
                print(f"  使用登录响应页面: {after_login_url}")
            except:
                pass
        
        # 如果还是没有，使用首页
        if not after_login_content:
            try:
                after_login_url = urljoin(base_url, '/')
                page_response = session.get(after_login_url, timeout=timeout)
                after_login_content = page_response.text
                print(f"  使用首页: {after_login_url}")
            except Exception as e:
                print(f"  ⚠ 无法获取登录后页面: {e}")
                after_login_content = login_response.text if hasattr(login_response, 'text') else ""
        
        # 验证登录是否成功（检查页面内容）
        login_verified = False
        if after_login_content:
            # 检查页面是否包含登录后的标识（而不是登录页面）
            login_success_indicators = [
                '学习空间', 'space', 'dashboard', '个人中心', 
                '我的课程', '退出登录', 'logout', '用户中心'
            ]
            login_page_indicators = ['登录', 'login', '用户名', 'password', '登录页面']
            
            content_lower = after_login_content.lower()
            has_success = any(indicator in content_lower for indicator in login_success_indicators)
            has_login_page = any(indicator in content_lower for indicator in login_page_indicators)
            
            # 如果URL是 /space/ 或包含成功标识，且不包含登录页面标识
            if '/space' in (after_login_url or '') or (has_success and not has_login_page):
                login_verified = True
        
        # 保存页面内容
        result['success'] = login_verified if after_login_content else login_api_success
        result['final_url'] = after_login_url if after_login_url else login_response.url
        result['content'] = after_login_content if after_login_content else (login_response.text if hasattr(login_response, 'text') else "")
        
        if result['content']:
            if login_response.encoding:
                result['encoding'] = login_response.encoding
            else:
                result['encoding'] = 'utf-8'
            
            if result['success']:
                print(f"\n✓ 登录成功并获取登录后页面！")
            else:
                print(f"\n⚠ 登录状态不确定，已获取页面内容")
            print(f"  最终URL: {result['final_url']}")
            print(f"  编码: {result['encoding']}")
            print(f"  页面大小: {len(result['content'].encode('utf-8'))} 字节")
            print(f"  文本长度: {len(result['content'])} 字符")
        
        # 保存HTML到文件
        if save_to_file and result['content']:
            filename = f"login_result_{int(time.time())}.html"
            try:
                with open(filename, 'w', encoding='utf-8') as f:
                    f.write(result['content'])
                print(f"\n✓ HTML内容已保存到文件: {filename}")
                result['saved_file'] = filename
            except Exception as e:
                print(f"\n⚠ 保存文件失败: {e}")
        
        # 显示页面内容（前1000个字符作为预览）
        if result['content']:
            print("\n" + "=" * 50)
            print("登录后页面内容预览（前1000字符）:")
            print("=" * 50)
            preview = result['content'][:1000]
            print(preview)
            if len(result['content']) > 1000:
                print(f"\n... (共 {len(result['content'])} 字符，完整内容已保存到文件)")
            print("=" * 50)
        
    except requests.exceptions.Timeout:
        result['error'] = "连接超时"
        print(f"✗ 连接失败: 超时（{timeout}秒）")
    except requests.exceptions.ConnectionError:
        result['error'] = "连接错误"
        print(f"✗ 连接失败: 无法连接到服务器")
    except requests.exceptions.RequestException as e:
        result['error'] = str(e)
        print(f"✗ 请求失败: {e}")
    except Exception as e:
        result['error'] = str(e)
        print(f"✗ 发生未知错误: {e}")
    finally:
        session.close()
    
    return result

def test_connection(url, timeout=10, show_content=True):
    """
    测试能否连通指定URL并获取页面内容
    
    参数:
        url: 要测试的URL
        timeout: 超时时间（秒）
        show_content: 是否显示页面内容（默认True）
    
    返回:
        dict: 包含连接状态和页面内容的字典
    """
    if not REQUESTS_AVAILABLE:
        result = {
            'url': url,
            'success': False,
            'status_code': None,
            'response_time': None,
            'content': None,
            'error': 'requests 模块未安装'
        }
        print(f"\n无法测试连接: {url}")
        print("-" * 50)
        print("✗ requests 模块未安装，无法进行网络测试")
        print("\n安装方法:")
        print("  1. 如果使用 conda 环境: conda activate base 然后 pip install requests")
        print("  2. 如果使用系统 Python: python3 -m pip install requests")
        print("  3. 或者直接使用 conda 环境中的 Python 运行此脚本")
        return result
    
    print(f"\n正在测试连接: {url}")
    print("-" * 50)
    
    result = {
        'url': url,
        'success': False,
        'status_code': None,
        'response_time': None,
        'content': None,
        'encoding': None,
        'error': None
    }
    
    try:
        start_time = time.time()
        # 设置请求头，模拟浏览器访问
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
        }
        response = requests.get(url, timeout=timeout, allow_redirects=True, headers=headers)
        end_time = time.time()
        
        result['success'] = True
        result['status_code'] = response.status_code
        result['response_time'] = round((end_time - start_time) * 1000, 2)  # 转换为毫秒
        result['final_url'] = response.url  # 最终URL（可能发生重定向）
        
        # 获取页面内容
        # 尝试自动检测编码
        if response.encoding:
            result['encoding'] = response.encoding
            try:
                result['content'] = response.text
            except UnicodeDecodeError:
                # 如果自动编码失败，尝试使用UTF-8
                result['content'] = response.content.decode('utf-8', errors='ignore')
                result['encoding'] = 'utf-8'
        else:
            # 如果没有指定编码，尝试UTF-8
            result['content'] = response.content.decode('utf-8', errors='ignore')
            result['encoding'] = 'utf-8'
        
        print(f"✓ 连接成功！")
        print(f"  状态码: {response.status_code}")
        print(f"  响应时间: {result['response_time']} ms")
        print(f"  最终URL: {result['final_url']}")
        print(f"  编码: {result['encoding']}")
        
        if response.status_code == 200:
            print(f"  页面大小: {len(response.content)} 字节")
            print(f"  文本长度: {len(result['content'])} 字符")
        
        # 显示页面内容
        if show_content and result['content']:
            print("\n" + "=" * 50)
            print("页面内容:")
            print("=" * 50)
            print(result['content'])
            print("=" * 50)
        
    except requests.exceptions.Timeout:
        result['error'] = "连接超时"
        print(f"✗ 连接失败: 超时（{timeout}秒）")
    except requests.exceptions.ConnectionError:
        result['error'] = "连接错误"
        print(f"✗ 连接失败: 无法连接到服务器")
    except requests.exceptions.RequestException as e:
        result['error'] = str(e)
        print(f"✗ 连接失败: {e}")
    except Exception as e:
        result['error'] = str(e)
        print(f"✗ 发生未知错误: {e}")
    
    return result

if __name__ == "__main__":
    # 登录信息
    login_url = "https://v.gaodun.com/login"
    username = "18643987165"
    password = "Maminjia95.,"
    
    # 已知的正确登录API端点（如果自动检测失败，可以手动指定）
    # login_api_url = "https://apigateway.gaodun.com/api/v1/ucenter/login"
    login_api_url = None  # None表示自动检测，或使用上面注释的URL
    
    print("=" * 50)
    print("网站登录脚本 - 模拟页面操作")
    print("=" * 50)
    print(f"操作流程:")
    print(f"  1. 访问登录页面: {login_url}")
    print(f"  2. 选择账户密码登录方式")
    print(f"  3. 填写用户名: {username}")
    print(f"  4. 填写密码并点击登录")
    print(f"  5. 获取登录后的页面内容")
    print("=" * 50)
    
    # 执行登录并获取登录后的页面内容
    login_result = login_and_get_content(login_url, username, password, timeout=15, save_to_file=True, api_url=login_api_url)
    
    print("\n" + "=" * 50)
    if login_result['success']:
        print("✓ 登录完成")
        if login_result['content']:
            print(f"✓ 已成功获取登录后的页面内容（{len(login_result['content'])} 字符）")
            if 'saved_file' in login_result:
                print(f"✓ HTML文件已保存: {login_result['saved_file']}")
    else:
        print("✗ 登录失败或状态不确定")
        if login_result['error']:
            print(f"  错误信息: {login_result['error']}")
        print("\n提示: 请检查:")
        print("  1. 用户名和密码是否正确")
        print("  2. 网络连接是否正常")
        print("  3. 网站登录页面结构是否发生变化")
