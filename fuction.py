import math


## 求Z值
def z(pc, tc, t, p):
    a1 = 0.31506237
    a2 = -1.0467099
    a3 = -0.57832729
    a4 = 0.53530771
    a5 = -0.61232032
    a6 = -0.10488813
    a7 = 0.68157001
    a8 = 0.68446549
    ppr = p / pc
    tpr = t / tc
    if 0.1 <= ppr <= 14:
        # Dranchuk-Purris-Robinson 方法
        luopr = 0.27 * ppr / tpr  # 初始值
        # 牛顿迭代30次
        for _ in range(30):
            # 计算当前函数值和导数
            fluopr = luopr - (0.27 * ppr) / tpr + (a1 + a2 / tpr + a3 / tpr ** 3) * luopr ** 2 + (
                    a4 + a5 / tpr) * luopr ** 3 + (a5 * a6 * luopr ** 6) / tpr + (a7 * luopr ** 3 / tpr ** 3) * (
                             1 + a8 * luopr ** 2) * math.exp(-a8 * luopr ** 2)
            dfluopr = 1 + (a1 + a2 / tpr + a3 / tpr ** 3) * 2 * luopr + (a4 + a5 / tpr) * 3 * luopr ** 2 + (
                    a5 * a6 / tpr) * (6 * luopr ** 5) + (a7 / tpr ** 3) * (
                              3 * luopr ** 2 + a8 * 3 * luopr ** 4) - a8 ** 2 * 2 * luopr ** 6 * math.exp(
                -a8 * luopr ** 2)
            luopr = luopr - fluopr / dfluopr
        z = 0.27 * ppr / (luopr * tpr)
    else:
        # Hall-Yarborough 方法
        tt = 1 / tpr
        yy = 0.06125 * ppr * tt * math.exp(-1.2 * (1 - tt) ** 2)  # 初始值
        # 牛顿迭代30次
        for _ in range(30):
            fhall = -0.06125 * ppr * tt * math.exp(-1.2 * (1 - tt) ** 2) + (yy + yy ** 2 + yy ** 3 - yy ** 4) / (
                    1 - yy) ** 3 - (14.76 * tt - 9.76 * tt ** 2 + 4.58 * tt ** 3) * yy ** 2 + (
                            90.7 * tt - 242.2 * tt ** 2 + 42.4 * tt ** 3) * yy ** (2.18 + 2.82 * tt)
            dhall = (1 + 4 * yy + 4 * yy ** 2 - 4 * yy ** 3 + yy ** 4) / (1 - yy) ** 4 - (
                    29.52 * tt - 19.52 * tt ** 2 + 9.16 * tt ** 3) * yy + (2.18 + 2.82 * tt) * (
                            90.7 * tt - 242.2 * tt ** 2 + 42.4 * tt ** 3) * yy ** (1.18 + 2.82 * tt)
            yy = yy - fhall / dhall
        z = 0.06125 * ppr * tt * math.exp(-1.2 * (1 - tt) ** 2) / yy
    return z


def bg(pc, tc, t, p):
    ## 天然气体积系数计算 + Dranchuk,Purris和Robinson法计算z
    zz = z(pc, tc, t, p)
    bg = 0.0003447 * zz * t / p
    return bg


def cg(pc, tc, t, p):
    # 天然气压缩系数
    # Dranchuk,Purvis和Robinson法计算z
    a1 = 0.31506237
    a2 = -1.0467099
    a3 = -0.57832729
    a4 = 0.53530771
    A5 = -0.61232032
    a6 = -0.10488813
    a7 = 0.68157001
    a8 = 0.68446549
    ppr = p / pc
    tpr = t / tc
    luopr = 0.27 * ppr / tpr
    for i in range(30, 0, -1):
        fluopr = luopr - (0.27 * ppr) / tpr + (a1 + a2 / tpr + a3 / tpr ** 3) * luopr ** 2 + (
                a4 + A5 / tpr) * luopr ** 3 + (A5 * a6 * luopr ** 6) / tpr + (a7 * luopr ** 3 / tpr ** 3) * (
                         1 + a8 * luopr ** 2) * math.exp(-a8 * luopr ** 2)
        dfluopr = 1 + (a1 + a2 / tpr + a3 / tpr ** 3) * (2 * luopr) + (a4 + A5 / tpr) * (3 * luopr ** 2) + (
                A5 * a6 / tpr) * (6 * luopr ** 5) + (a7 / tpr ** 3) * (
                          3 * luopr ** 2 + a8 * (3 * luopr ** 4)) - a8 ** 2 * (2 * luopr ** 6) * math.exp(
            -a8 * luopr ** 2)
        luopr = luopr - fluopr / dfluopr
    z = 0.27 * ppr / (luopr * tpr)
    # 杨继盛“采气工艺基础”（旧）34页
    dzlt = (a1 + a2 / tpr + a3 / tpr ** 3) + 2 * (a4 + A5 / tpr) * luopr + 5 * A5 * a6 * luopr ** 4 / tpr + (
            2 * a7 * luopr / tpr ** 3) * (1 + a8 * luopr ** 2 - a8 ** 2 * luopr ** 4) * math.exp(-a8 * luopr ** 2)
    cpr = 1 / ppr - (0.27 / (z ** 2 * tpr)) * (dzlt / (1 + luopr * dzlt / z))
    cg = cpr / pc
    return cg


def llupr(rg, pc, tc, t, p):
    # Dranchuk,Purris和Robinson法计算z
    # Lee,Gonzalez和Eakin法计算密度ρ
    # 杨继盛“采气工艺基础”（旧）40页
    # Dranchuk,Purris和Robinson法计算z
    zz = z(pc, tc, t, p)
    lluopr = 3.4844 * p * rg / (zz * t)  # lluopr----密度
    return lluopr

def zhandu(rg, pc, tc, t, p, yn2, yco2, yh2s):
    # Lee,Gonzalez和Eakin法计算粘度μ,单位cp
    # 杨继盛“采气工艺基础”（旧）40页
    # 酸性气体修正（1986年）

    kn2 = yn2 * (0.00005 * rg + 0.000047) * 100
    kco2 = yco2 * (0.000078 * rg + 0.00001) * 100
    kh2s = yh2s * (0.000058 * rg - 0.000018) * 100

    k = (0.0001 * (9.4 + 0.02 * 28.97 * rg) * (9 * t / 5) ** 1.5) / (
            209 + 19 * 28.97 * rg + 9 * t / 5) + kn2 + kco2 + kh2s
    x = 3.5 + 986 / (9 * t / 5) + 0.01 * 28.97 * rg
    y = 2.4 - 0.2 * x

    lllupr = llupr(rg, pc, tc, t, p)  # 计算密度
    zhandu = k * math.exp(x * (lllupr ** y))  # zhandu----粘度
    return zhandu
#####################需要核实！！！！！！！！！！！！！！！
def zp(pc, tc, t, p):
    # z()函数为Dranchuk,Purris和Robinson法计算z
    for i in range(30, 0, -1):
        zz = z(pc, tc, t, p)
        p = p * zz
        zp = p
    return  zp

def ft(rg, pc, tc, tts, tws, p, yn2, yco2, yh2s, d1, d3, q, ee):
    # Lee,Gonzalez和Eakin法计算粘度μ
    # 杨继盛“采气工艺基础”（旧）40页
    # 酸性气体修正（1986年）
    # 套管生产
    kn2 = yn2 * (0.00005 * rg + 0.000047) * 100
    kco2 = yco2 * (0.000078 * rg + 0.00001) * 100
    kh2s = yh2s * (0.000058 * rg - 0.000018) * 100
    t = (tts + tws) / 2
    k = (0.0001 * (9.4 + 0.02 * 28.97 * rg) * (9 * t / 5) ^ 1.5) / (
                209 + 19 * 28.97 * rg + 9 * t / 5) + kn2 + kco2 + kh2s
    x = 3.5 + 986 / (9 * t / 5) + 0.01 * 28.97 * rg
    y = 2.4 - 0.2 * x

    #Dranchuk,Purris和Robinson法计算z
    zz = z(pc, tc, t, p)
    # lluopr----密度
    lluopr = 0.0014926 * (144.9275 * p) * (28.97 * rg) / (zz * (9 * t / 5))
    # zhandu----粘度
    zhandu = k * math.exp(x * lluopr ** y)
    # Jain法计算摩阻系数f
    # 杨继盛“采气工艺基础”（旧）3-30页
    re = 179.39789 * q * rg / ((d1 + d3) * zhandu)
    # d1----油管直径
    ft = 1 / (1.14 - 2 * math.log10(ee / (d1 + d3) + 21.25 / re  ** 0.9)/math.log10(10)) ** 2
    return ft

def fy(rg, pc, tc, tts, tws, p, yn2, yco2, yh2s, d1, q, ee):
    # Lee, Gonzalez和Eakin法计算粘度μ
    # 杨继盛“采气工艺基础”（旧）40页
    # 酸性气体修正（1986年）
    # 油管生产
    # ee----粗糙系数

    kn2 = yn2 * (0.00005 * rg + 0.000047) * 100
    kco2 = yco2 * (0.000078 * rg + 0.00001) * 100
    kh2s = yh2s * (0.000058 * rg - 0.000018) * 100

    t = (tts + tws) / 2
    k = (0.0001 * (9.4 + 0.02 * 28.97 * rg) * (9 * t / 5) ** 1.5) / (209 + 19 * 28.97 * rg + 9 * t / 5) + kn2 + kco2 + kh2s
    x = 3.5 + 986 / (9 * t / 5) + 0.01 * 28.97 * rg
    y = 2.4 - 0.2 * x

    # Dranchuk, Purris和Robinson法计算z
    zz = z(pc, tc, t, p)  # 假设z函数已经定义好

    lluopr = 0.0014926 * (144.9275 * p) * (28.97 * rg) / (zz * (9 * t / 5))  # lluopr----密度
    zhandu = k * math.exp(x * (lluopr ** y))  # zhandu----粘度
    # Jain法计算摩阻系数f
    # 杨继盛“采气工艺基础”（旧）3-30页
    re = 179.39789 * q * rg / (d1 * zhandu)  # d1----油管直径
    fy = 1 / (1.14 - 2 * math.log(ee / d1 + 21.25 / (re ** 0.9)) / math.log(10)) ** 2
    return fy

def pws(rg, pc, tc, h, tts, tws, pts):
    # 平均温度和平均压缩系数计算法计算井底压力(静气柱）

    pws = pts + pts * h / 12192

    for i in range(30, 0, -1):
        p = (pts + pws_value) / 2
        t = (tts + tws) / 2

        # Dranchuk, Purris和Robinson法计算z
        zz = z(pc, tc, t, p)  # 假设z函数已经定义好
        pws = pts * math.exp(0.03415 * rg * h / (zz * t))

    return pws

def pts(rg, pc, tc, h, tts, tws, pws):
    # 平均温度和平均压缩系数计算法计算井口压力(静气柱）

    pts = pws - pws * h / 12192

    for i in range(100, 0, -1):
        p = (pts + pws) / 2
        t = (tts + tws) / 2

        # Dranchuk, Purris和Robinson法计算z
        zz = z(pc, tc, t, p)  # 假设z函数已经定义好

        pts_value = pws / math.exp(0.03415 * rg * h / (zz * t))

    return pts

def pwf(rg, pc, tc, yn2, yco2, yh2s, d1, d2, d3, h1, h2, h3, h, tts, tws, q, ptf, ee):
    # 平均温度和平均压缩系数计算油管采气时井底流动压力Pwf

    # d1、h1---第1段油管直径和下入长度
    # d2、h2---第2段油管直径和下入长度
    # d3、h3---产层直径和油管底部至中部井深的长度

    tj = 0.8742 * q + 20.22 + 273.15  # tj----川东井口温度计算经验公式

    # 计算第1段压力
    t1 = h1 * (tws - tj) / h + tj  # t1----第1段油管底部处温度
    t = (tj + t1) / 2

    if d1 > 0:
        # 赋初值
        pwf1 = ptf + ptf * h1 / 12192  # pwf1----第1段油管底部处压力

        for i1 in range(30, 0, -1):
            p = (pwf1 + ptf) / 2

            # Dranchuk, Purris和Robinson法计算z
            zz = z(pc, tc, t, p)  # 假设z函数已经定义好

            # Jain法计算摩阻系数f
            # 杨继盛“采气工艺基础”（旧）3-30页
            ffy = fy(rg, pc, tc, tts, tws, p, yn2, yco2, yh2s, d1, q, ee)
            S = 0.03415 * rg * h1 / (t * zz)
            pwf1 = math.sqrt(ptf ** 2 * math.exp(2 * S) + (1.324 * 0.0000000001 * ffy * (q * t * zz) ** 2 * (math.exp(2 * S) - 1)) / d1 ** 5)
    else:
        pwf1 = ptf

    # 计算第2段压力
    t2 = (h1 + h2) * (tws - tj) / h + tj
    t = (t1 + t2) / 2

    if d2 > 0:
        # 赋初值
        pwf2 = pwf1 + pwf1 * h2 / 12192  # pwf2----第2段油管底部处压力

        for i2 in range(30, 0, -1):
            p = (pwf1 + pwf2) / 2

            # Dranchuk, Purris和Robinson法计算z
            zz = z(pc, tc, t, p)  # 假设z函数已经定义好

            # Jain法计算摩阻系数f
            # 杨继盛“采气工艺基础”（旧）3-30页
            d1 = d2
            ffy = fy(rg, pc, tc, tts, tws, p, yn2, yco2, yh2s, d1, q, ee)

            S = 0.03415 * rg * h2 / (t * zz)
            pwf2 = math.sqrt(pwf1 ** 2 * math.exp(2 * S) + (1.324 * 0.0000000001 * ffy * (q * t * zz) ** 2 * (math.exp(2 * S) - 1)) / d1 ** 5)
    else:
        pwf2 = pwf1

    # 计算井底压力
    t2 = (h1 + h2) * (tws - tj) / h + tj
    t = (t2 + tws) / 2

    if d3 > 0:
        # 赋初值
        pwf = pwf2 + pwf2 * h3 / 12192

        for i3 in range(30, 0, -1):
            p = (pwf2 + pwf) / 2

            # Dranchuk, Purris和Robinson法计算z
            zz = z(pc, tc, t, p)  # 假设z函数已经定义好

            # Jain法计算摩阻系数f
            # 杨继盛“采气工艺基础”（旧）3-30页
            d1 = d3
            ffy = fy(rg, pc, tc, tts, tws, p, yn2, yco2, yh2s, d1, q, ee)

            S = 0.03415 * rg * (h - h1 - h2) / (t * zz)
            pwf = math.sqrt(pwf2 ** 2 * math.exp(2 * S) + (1.324 * 0.0000000001 * ffy * (q * t * zz) ** 2 * (math.exp(2 * S) - 1)) / d1 ** 5)
    else:
        pwf = pwf2

    return pwf

def ptf(rg, pc, tc, yn2, yco2, yh2s, d4, d2, d3, h1, h2, h3, h, tts, tws, q, pwf, ee):
    # 平均温度和平均压缩系数计算油管采气时井口流动压力Pwf

    tj = 0.8742 * q + 20.22 + 273.15  # tj----川东井口温度计算经验公式

    # 计算第3段压力（从井底向上算）
    t2 = (h1 + h2) * (tws - tj) / h + tj  # t2----油管底部处温度
    t = (tws + t2) / 2

    if d3 > 0:
        # 赋初值
        ptf2 = pwf - pwf * (h - h1 - h2) / 121920

        for i3 in range(30, 0, -1):
            p = (pwf + ptf2) / 2

            # Dranchuk, Purris和Robinson法计算z
            zz = z(pc, tc, t, p)  # 假设z函数已经定义好

            # Jain法计算摩阻系数f
            # 杨继盛“采气工艺基础”（旧）3-30页
            d1 = d3
            ffy = fy(rg, pc, tc, tts, tws, p, yn2, yco2, yh2s, d1, q, ee)

            S = 0.03415 * rg * (h - h1 - h2) / (t * zz)
            ptf2 = math.sqrt((pwf ** 2 - (1.324 * 0.0000000001 * ffy * (q * t * zz) ** 2 * (math.exp(2 * S) - 1)) / d1 ** 5) / math.exp(2 * S))
    else:
        ptf2 = pwf

    # 计算第2段压力（从井底向上算）
    t1 = h1 * (tws - tj) / h + tj  # t1----第1段油管连接处温度
    t = (t1 + t2) / 2

    if d2 > 0:
        # 赋初值
        ptf1 = ptf2 - ptf2 * h2 / 121920

        for i2 in range(30, 0, -1):
            p = (ptf1 + ptf2) / 2

            # Dranchuk, Purris和Robinson法计算z
            zz = z(pc, tc, t, p)  # 假设z函数已经定义好

            # Jain法计算摩阻系数f
            # 杨继盛“采气工艺基础”（旧）3-30页
            d1 = d2
            ffy = fy(rg, pc, tc, tts, tws, p, yn2, yco2, yh2s, d1, q, ee)

            S = 0.03415 * rg * h2 / (t * zz)
            ptf1 = math.sqrt((ptf2 ** 2 - (1.324 * 0.0000000001 * ffy * (q * t * zz) ** 2 * (math.exp(2 * S) - 1)) / d1 ** 5) / math.exp(2 * S))
    else:
        ptf1 = ptf2

    # 计算井口压力（从井底向上算）
    t = (t1 + tj) / 2

    if d4 > 0:
        # 赋初值
        ptf_value = ptf1 - ptf1 * h1 / 121920

        for i1 in range(30, 0, -1):
            p = (ptf1 + ptf_value) / 2

            # Dranchuk, Purris和Robinson法计算z
            zz = z(pc, tc, t, p)  # 假设z函数已经定义好

            # Jain法计算摩阻系数f
            # 杨继盛“采气工艺基础”（旧）3-30页
            d1 = d4
            ffy = fy(rg, pc, tc, tts, tws, p, yn2, yco2, yh2s, d1, q, ee)

            S = 0.03415 * rg * h1 / (t * zz)
            ptf_value = math.sqrt((ptf1 ** 2 - (1.324 * 0.0000000001 * ffy * (q * t * zz) ** 2 * (math.exp(2 * S) - 1)) / d1 ** 5) / math.exp(2 * S))
    else:
        ptf_value = ptf1

    return ptf_value

def pwfh(rg, pc, tc, yn2, yco2, yh2s, d1, d2, d3, h1, h2, h3, h, tts, tws, q, ptf, ee):
    # 平均温度和平均压缩系数计算环空采气时井底流动压力Pwfh

    tj = 0.8742 * q + 20.22 + 273.15  # tj----川东井口温度计算经验公式

    # 计算第1段压力
    t1 = h1 * (tws - tj) / h + tj  # t1----第1段油管底部处温度
    t = (tj + t1) / 2

    # 赋初值
    pwfh1 = ptf + ptf * h / 12192

    for i1 in range(1, 0, -1):
        p = (pwfh1 + ptf) / 2

        # Dranchuk, Purris和Robinson法计算z
        zz = z(pc, tc, t, p)  # 假设z函数已经定义好

        # Jain法计算摩阻系数f
        # 杨继盛“采气工艺基础”（旧）3-30页
        fft = ft(rg, pc, tc, tts, tws, p, yn2, yco2, yh2s, d1, d3, q, ee)
        S = 0.03415 * rg * h1 / (t * zz)
        pwfh1 = math.sqrt(ptf ** 2 * math.exp(2 * S) + (1.324 * 0.0000000001 * fft * (q * t * zz) ** 2 * (math.exp(2 * S) - 1)) / ((d3 - d1) ** 3 * (d3 + d1) ** 2))

    # 计算第2段压力
    t2 = (h1 + h2) * (tws - tj) / h + tj
    t = (t1 + t2) / 2

    # 赋初值
    pwfh2 = pwfh1 + pwfh1 * h / 12192

    for i2 in range(1, 0, -1):
        p = (pwfh1 + pwfh2) / 2

        # Dranchuk, Purris和Robinson法计算z
        zz = z(pc, tc, t, p)  # 假设z函数已经定义好

        # Jain法计算摩阻系数f
        # 杨继盛“采气工艺基础”（旧）3-30页
        d1 = d2
        fft = ft(rg, pc, tc, tts, tws, p, yn2, yco2, yh2s, d1, d3, q, ee)
        S = 0.03415 * rg * h2 / (t * zz)
        pwfh2 = math.sqrt(pwfh1 ** 2 * math.exp(2 * S) + (1.324 * 0.0000000001 * fft * (q * t * zz) ** 2 * (math.exp(2 * S) - 1)) / ((d3 - d1) ** 3 * (d3 + d1) ** 2))

    # 计算井底压力
    t2 = (h1 + h2) * (tws - tj) / h + tj
    t = (t2 + tws) / 2

    # 赋初值
    pwfh = pwfh2 + pwfh2 * h3 / 12192

    for i3 in range(1, 0, -1):
        p = (pwfh2 + pwfh) / 2

        # Dranchuk, Purris和Robinson法计算z
        zz = z(pc, tc, t, p)  # 假设z函数已经定义好

        # Jain法计算摩阻系数f
        # 杨继盛“采气工艺基础”（旧）3-30页
        d1 = d3
        ffy = fy(rg, pc, tc, tts, tws, p, yn2, yco2, yh2s, d1, q, ee)
        S = 0.03415 * rg * (h - h1 - h2) / (t * zz)
        pwfh = math.sqrt(pwfh2 ** 2 * math.exp(2 * S) + (1.324 * 0.0000000001 * ffy * (q * t * zz) ** 2 * (math.exp(2 * S) - 1)) / d1 ** 5)

    return pwfh

def ptfh(rg, pc, tc, yn2, yco2, yh2s, d4, d2, d3, h1, h2, h3, h, tts, tws, q, pwfh, ee):
    # 平均温度和平均压缩系数计算环空采气时井口流动压力Pwfh

    tj = 0.8742 * q + 20.22 + 273.15  # tj----川东井口温度计算经验公式

    # 计算第3段压力（从井底向上算）
    t2 = (h1 + h2) * (tws - tj) / h + tj  # t2----油管底部处温度
    t = (tws + t2) / 2

    # 赋初值
    ptfh2 = pwfh - pwfh * (h - h1 - h2) / 121920

    for i3 in range(30, 0, -1):
        p = (pwfh + ptfh2) / 2

        # Dranchuk, Purris和Robinson法计算z
        zz = z(pc, tc, t, p)  # 假设z函数已经定义好

        # Jain法计算摩阻系数f
        # 杨继盛“采气工艺基础”（旧）3-30页
        d1 = d3
        ffy = fy(rg, pc, tc, tts, tws, p, yn2, yco2, yh2s, d1, q, ee)

        S = 0.03415 * rg * (h - h1 - h2) / (t * zz)
        ptfh2 = math.sqrt((pwfh ** 2 - (1.324 * 0.0000000001 * ffy * (q * t * zz) ** 2 * (math.exp(2 * S) - 1)) / d1 ** 5) / math.exp(2 * S))

    # 计算第2段压力（从井底向上算）
    t1 = h1 * (tws - tj) / h + tj  # t1----第1段油管连接处温度
    t = (t1 + t2) / 2

    # 赋初值
    ptfh1 = ptfh2 - ptfh2 * h2 / 121920

    for i2 in range(30, 0, -1):
        p = (ptfh1 + ptfh2) / 2

        # Dranchuk, Purris和Robinson法计算z
        zz = z(pc, tc, t, p)  # 假设z函数已经定义好

        # Jain法计算摩阻系数f
        # 杨继盛“采气工艺基础”（旧）3-30页
        d1 = d2
        fft = ft(rg, pc, tc, tts, tws, p, yn2, yco2, yh2s, d1, d3, q, ee)

        S = 0.03415 * rg * h2 / (t * zz)
        ptfh1 = math.sqrt((ptfh2 ** 2 - (1.324 * 0.0000000001 * fft * (q * t * zz) ** 2 * (math.exp(2 * S) - 1)) / ((d3 - d2) ** 3 * (d2 + d3) ** 2)) / math.exp(2 * S))

    # 计算井口压力
    t = (t1 + tj) / 2

    # 赋初值
    ptfh_value = ptfh1 - ptfh1 * h1 / 121920

    for i3 in range(30, 0, -1):
        p = (ptfh_value + ptfh1) / 2

        # Dranchuk, Purris和Robinson法计算z
        zz = z(pc, tc, t, p)  # 假设z函数已经定义好

        # Jain法计算摩阻系数f
        # 杨继盛“采气工艺基础”（旧）3-30页
        d1 = d4
        fft = ft(rg, pc, tc, tts, tws, p, yn2, yco2, yh2s, d1, d3, q, ee)
        S = 0.03415 * rg * h1 / (t * zz)
        ptfh_value = math.sqrt((ptfh1 ** 2 - (1.324 * 0.0000000001 * fft * (q * t * zz) ** 2 * (math.exp(2 * S) - 1)) / ((d3 - d1) ** 3 * (d1 + d3) ** 2)) / math.exp(2 * S))

    return ptfh_value

def qkp(rg, pc, tc, t, p, d):
    # 气井连续带液临界流量计算，采气工程，杨川东，P102

    # qkp----气井连续带液临界流量,104m3/d
    # t----井底温度，K
    # p----井底流动压力，MPa
    # d----油管内径，cm

    zz = z(pc, tc, t, p)  # 假设z函数已经定义好
    d = d * 100  # 将油管内径从米转换为厘米
    qkp_value = 0.0648 * (rg * zz * t) ** (-0.5) * (10553 - (34158 * rg * p) / (zz * t)) ** 0.25 * p ** 0.5 * d ** 2

    return qkp_value

def di(rg, pc, tc, t, p, q):
    # 气井连续带液合理油管直径计算，采气工程，杨川东，P102

    # q----气井连续带液临界流量,104m3/d
    # t----井底温度，K
    # p----井底流动压力，MPa
    # di----油管内径，m

    zz = z(pc, tc, t, p)  # 假设z函数已经定义好
    di_value = 0.01 * 1.2423 * (rg * zz * t) ** 0.25 * (10553 - (34158 * rg * p) / (zz * t)) ** (-1 / 8) * p ** (-0.25) * q ** 0.5

    return di_value

rlt = zp(pc=4.69, tc=194.4, t=376.04, p=25.68)
print(rlt)