#import library
import pandas as pd
import numpy as np
import math
import datetime

############################################################################################################################
#定数
############################################################################################################################
Air_Cp   = 1.005                                                                        #空気の定圧比熱　   kJ / (kg・K)
Vap_Cp   = 1.846                                                                        #水蒸気の定圧比熱   KJ / (kg・K)
Vap_L    = 2501.1                                                                       #水蒸気の蒸発潜熱   KJ / (kg・K)

Sigma    = 4.88e-8                                                                      #ステファン・ボルツマン定数
Solar_I = 1365                                                                          #太陽定数           W / m2

capa_air = lambda v: v * 1.005 * 1.2 * 1000                                             #容積v[m3]の空気の熱容量　J/K

############################################################################################################################
#湿り空気の状態
############################################################################################################################
get_rho = lambda sita:   353.25 / (sita + 273.15)                                       #空気の密度rhoを返す  

T      = lambda t: t + 273.15                                                           #絶対温度

T_dash = lambda t: T(100.0) / T(t)                                                      #飽和水蒸気の計算に用いる温度

Wh_to_MJ   = lambda v: v * 3.6 / 1000                                                   #Wh->MJ
MJ_to_Wh   = lambda v: v * 1000 / 3.6

f  = lambda t:    - 7.90295 * (T_dash(t) - 1) \
                  + 5.02808 * np.log10(T_dash(t)) \
                  - 1.3816e-7 * (np.power(10, 11.344  * (1 - T_dash(t))) - 1) \
                  + 8.1328e-3 * (np.power(10, -3.4919 * (T_dash(t- 1))) - 1) \
                  + np.log10(1013.246)

f_from_x = lambda x: (x / 1000 * 101325) / (x / 1000 + 0.622)

ps = lambda t:    np.power(10, f(t)) * 100                                                      #飽和水蒸気圧

e  = lambda t, h: h / 100 * ps(t)                                                               #水蒸気圧

x  = lambda t, h: 0.622 * (e(t, h) / (101325 - e(t, h)))                                        #絶対湿度

hs = lambda t:    Air_Cp * t                                                                    #エンタルピ（顕熱）

hl = lambda t, h: x(t, h) * (Vap_L + Vap_Cp * t)                                                #エンタルピ（潜熱）

ht = lambda t, h: hs(t) + hl(t, h)                                                              #エンタルピ（全熱）

############################################################################################################################
#風向
############################################################################################################################
def make_wind(d, s, c_in = 0.7, c_out = -0.55, c_horizontal = -0.90):
    df = pd.DataFrame(index = d.index) 

    df['w_d']    = d
    df['w_s']    = s
    df['風向']   = d * 22.5

    df['風速_E'] = np.sin(np.radians(d * 22.5)) * s
    df['風速_S'] = - np.cos(np.radians(d * 22.5)) * s
    df['風速_W'] = - np.sin(np.radians(d * 22.5)) * s
    df['風速_N'] = np.cos(np.radians(d * 22.5)) * s

    df.loc[df['風速_E'] >= 0, '風圧_E'] =  1.2 / 2 * c_in  * df['風速_E'] ** 2
    df.loc[df['風速_E'] < 0,  '風圧_E'] = -1.2 / 2 * c_out * df['風速_E'] ** 2
    df.loc[df['風速_S'] >= 0, '風圧_S'] =  1.2 / 2 * c_in  * df['風速_S'] ** 2
    df.loc[df['風速_S'] < 0,  '風圧_S'] = -1.2 / 2 * c_out * df['風速_S'] ** 2
    df.loc[df['風速_W'] >= 0, '風圧_W'] =  1.2 / 2 * c_in  * df['風速_W'] ** 2
    df.loc[df['風速_W'] < 0,  '風圧_W'] = -1.2 / 2 * c_out * df['風速_W'] ** 2
    df.loc[df['風速_N'] >= 0, '風圧_N'] =  1.2 / 2 * c_in  * df['風速_N'] ** 2
    df.loc[df['風速_N'] < 0,  '風圧_N'] = -1.2 / 2 * c_out * df['風速_N'] ** 2
    df['風速_H']                        =  1.2 / 2 * c_horizontal * df['w_s']

    #0:無風,  1:NNE,  2:NE,   3:ENE, 4:E,    5:ESE,  6:SE,   7:SSE,  8:S
                                          #9:SSW,  10:SW,  11:WSW, 12:W,  13:WNW, 14:NW,  15:NNW, 16:N
    wind_pressure = {
        '風圧_E': df['風圧_E'],
        '風圧_S': df['風圧_S'],  
        '風圧_W': df['風圧_W'],
        '風圧_N': df['風圧_N'],
        '風圧_H': df['風圧_H']
    }
    return df, wind_pressure

############################################################################################################################
#夜間放射
############################################################################################################################
rn = lambda t, h: (94.21 + 39.06 * np.sqrt(e(t, h) / 100) \
                   - 0.85 * Sigma * np.power(T(t), 4)) * 4.187 / 1000                           #夜間放射 MJ/m2

def make_nocturnal(**kwargs):
    
    if('t' in kwargs): 
        t = kwargs['t']
        h = kwargs['h']
        df = pd.DataFrame(index = t.index) 
        df['nocturnal'] = MJ_to_Wh(rn(t, h))
    elif('n_r' in kwargs):
        n_r = kwargs['n_r']
        df = pd.DataFrame(index = n_r.index)
        df['nocturnal'] = n_r
    else:                   
        raise Exception('ERROR: 温度 t がありません。夜間放射量 n_r もありません。')

    nocturnal = {
        '夜間_H':       df['nocturnal'],
        '夜間_V':       df['nocturnal'] * 0.5
    }

    return df, nocturnal

############################################################################################################################
#日射量
############################################################################################################################

#直散分離 erbs法
Kt = lambda IG, alt:     IG / (Wh_to_MJ(Solar_I) * np.sin(np.radians(alt)))                               #晴天指数

def Id(IG, kt):                                                                                 #水平面拡散日射量
    s_Id = np.zeros(len(kt))
    for i, k in enumerate(kt):
        if   k <= 0.22:                 s_Id[i] = IG[i] * (1 - 0.09 * k)
        elif (0.22 < k) & (k <= 0.80):  s_Id[i] = IG[i] * (0.9511 -  0.1604 * k \
                                                                  +  4.388  * np.power(k, 2) \
                                                                  - 16.638  * np.power(k, 3) \
                                                                  + 12.336  * np.power(k, 4))
        elif 0.80 < k:                  s_Id[i] = 0.365 * IG[i] 
    return s_Id

def Ib(IG, Id, alt):                                                                            #法線面直達日射量
    s_Ib = np.zeros(len(Id))
    for i, id in enumerate(Id):
        s_Ib[i] =  (IG[i] - Id[i]) / np.sin(np.radians(alt[i]))
        if (alt[i] < 10.0) & (s_Ib[i] > IG[i]):  s_Ib[i] = IG[i]
    return s_Ib

#太陽位置の計算
delta_d = lambda N: (180 / np.pi) * (0.006322 \
                                     - 0.405748 * np.cos(2 * np.pi * N / 366 + 0.153231) \
                                     - 0.005880 * np.cos(4 * np.pi * N / 366 + 0.207099) \
                                     - 0.003233 * np.cos(6 * np.pi * N / 366 + 0.620129))                               #太陽の赤緯

e_d     = lambda N: -0.000279 + 0.122772 * np.cos(2 * np.pi * N / 366 + 1.498311) \
                              - 0.165458 * np.cos(4 * np.pi * N / 366 - 1.261546) \
                              - 0.005354 * np.cos(6 * np.pi * N / 366 - 1.1571)                                         #太陽の均時差

T_d_t   = lambda H, ed, L       : (H  + ed - 12.0) * 15.0 + (L - 135.0)                                                 #太陽の時角

sin     = lambda v              : np.sin(np.radians(v))
cos     = lambda v              : np.cos(np.radians(v))

sin_hs  = lambda L, dd, tdt     : sin(L) * sin(dd) + cos(L) * cos(dd) * cos(tdt)                                        #仰角の正弦
sin_AZs = lambda dd, tdt, c_h   : cos(dd) * sin(tdt) / c_h                                                              #方位角の正弦
cos_AZs = lambda s_h, L, dd, c_h: (s_h * sin(L) - sin(dd)) / (c_h * cos(L))                                             #方位角の余弦

def sun_loc(idx, lat = 36.00, lon = 140.00, td = -0.5):
    df = pd.DataFrame(index = idx)
    df['N']       = [(i - datetime.datetime(i.year, 1, 1)).days + 1.5 for i in idx]                                     #元日からの通し日数
    df['H']       = idx.strftime("%H").astype('float64') + idx.strftime("%M").astype('float64')  / 60 + td              #時刻
    df['delta_d'] = delta_d(df['N'])                                                                                    #太陽の赤緯
    df['e_d']     = e_d(df['N'])                                                                                        #太陽の均時差
    df['T_d_t']   = T_d_t(df['H'], df['e_d'], lon)                                                                      #太陽の時角
    df['sin_hs']  = sin_hs(lat, df['delta_d'], df['T_d_t'])                                                             #太陽高度の正弦　sin
    df['cos_hs']  = np.sqrt(1 - np.power(df['sin_hs'], 2))                                                              #太陽高度の余弦　cos
    df['hs']      = np.degrees(np.arcsin(df['sin_hs']))                                                                 #太陽高度　°

    df['sin_AZs'] = sin_AZs(df['delta_d'], df['T_d_t'], df['cos_hs'])                                                   #太陽方位角の正弦　sin
    df['cos_AZs'] = cos_AZs(df['sin_hs'], lat, df['delta_d'], df['cos_hs'])                                             #太陽方位角の余弦　cos

    df.loc[(df['sin_AZs'] <   0) & (df['cos_AZs'] >  0), 'AZs'] = np.degrees(np.arctan(df['sin_AZs'] / df['cos_AZs']))
    df.loc[(df['sin_AZs'] >   0) & (df['cos_AZs'] >  0), 'AZs'] = np.degrees(np.arctan(df['sin_AZs'] / df['cos_AZs']))
    df.loc[(df['sin_AZs'] >   0) & (df['cos_AZs'] <  0), 'AZs'] = np.degrees(np.arctan(df['sin_AZs'] / df['cos_AZs'])) + 180
    df.loc[(df['sin_AZs'] <   0) & (df['cos_AZs'] <  0), 'AZs'] = np.degrees(np.arctan(df['sin_AZs'] / df['cos_AZs'])) - 180
    df.loc[(df['sin_AZs'] ==  1) & (df['cos_AZs'] == 0), 'AZs'] = 90
    df.loc[(df['sin_AZs'] == -1) & (df['cos_AZs'] == 0), 'AZs'] = -90                                                   #太陽方位角　°

    return df

def astro_sun_loc(idx, lat = '36 00 00.00', lon = '140 00 00.00', td = -0.5):
    import astropy.time
    import astropy.units as u
    from astropy.coordinates import get_sun
    from astropy.coordinates import AltAz
    from astropy.coordinates import EarthLocation

    loc = EarthLocation(lat = lat, lon = lon)
    time = astropy.time.Time(idx) + (-9 + td ) * u.hour
    sun = get_sun(time).transform_to(AltAz(obstime = time, location = loc))

    df = pd.DataFrame(index = idx)

    df['s_alt'] = np.array([np.sin(s.alt) for s in sun]).astype('float64')
    df['c_alt'] = np.array([np.cos(s.alt) for s in sun]).astype('float64')
    df['alt']   = np.degrees(np.arcsin(df['s_alt']))

    df['s_az']  = np.array([np.sin(s.az) for s in sun]).astype('float64')
    df['c_az']  = np.array([np.cos(s.az) for s in sun]).astype('float64') 
    df.loc[(df['s_az'] < 0) & (df['c_az'] > 0), 'az'] = np.degrees(np.arctan(df['s_az'] / df['c_az'])) + 180
    df.loc[(df['s_az'] > 0) & (df['c_az'] > 0), 'az'] = np.degrees(np.arctan(df['s_az'] / df['c_az'])) - 180
    df.loc[(df['s_az'] > 0) & (df['c_az'] < 0), 'az'] = np.degrees(np.arctan(df['s_az'] / df['c_az']))
    df.loc[(df['s_az'] < 0) & (df['c_az'] < 0), 'az'] = np.degrees(np.arctan(df['s_az'] / df['c_az']))

    return df

def sep_direct_diffuse(s_ig, s_hs):
    df_i = pd.concat([s_ig, s_hs], axis = 1)
    df_i.columns = ['IG', 'hs']
    
    df_i['Kt'] = Kt(Wh_to_MJ(df_i['IG']), df_i['hs'])
    df_i['Id'] = MJ_to_Wh(Id(Wh_to_MJ(df_i['IG']), df_i['Kt']))
    df_i['Ib'] = MJ_to_Wh(Ib(Wh_to_MJ(df_i['IG']), Wh_to_MJ(df_i['Id']), df_i['hs'])) 
    return df_i[['Kt', 'Id', 'Ib']]

def direc_solar(s_ib, s_id, s_sin_hs, s_cos_hs, s_hs, s_sin_AZs, s_cos_AZs, s_AZs):
    df_i = pd.concat([s_ib, s_id, s_sin_hs, s_cos_hs, s_hs, s_sin_AZs, s_cos_AZs, s_AZs], axis = 1)
    df_i.columns = ['Ib', 'Id', 'sin_hs', 'cos_hs', 'hs', 'sin_AZs', 'cos_AZs', 'AZs']

    df_i.loc[(df_i['hs'] > 0) & (-180 < df_i['AZs']) & (df_i['AZs'] < 0),   'Ib_E'] = -1 * df_i['Ib'] * df_i['cos_hs'] * df_i['sin_AZs']    #東面   E
    df_i.loc[(df_i['hs'] > 0) & (-90  < df_i['AZs']) & (df_i['AZs'] < 90),  'Ib_S'] =      df_i['Ib'] * df_i['cos_hs'] * df_i['cos_AZs']    #南面   S
    df_i.loc[(df_i['hs'] > 0) & (0    < df_i['AZs']) & (df_i['AZs'] < 180), 'Ib_W'] =      df_i['Ib'] * df_i['cos_hs'] * df_i['sin_AZs']    #西面   W
    df_i.loc[(df_i['hs'] > 0) & (-180 < df_i['AZs']) & (df_i['AZs'] < -90), 'Ib_N'] = -1 * df_i['Ib'] * df_i['cos_hs'] * df_i['cos_AZs']    #北面   N
    df_i.loc[(df_i['hs'] > 0) & (  90 < df_i['AZs']) & (df_i['AZs'] < 180), 'Ib_N'] = -1 * df_i['Ib'] * df_i['cos_hs'] * df_i['cos_AZs']    #北面   N
    df_i.loc[df_i['hs'] > 0,  'Ib_H'] = df_i['Ib'] * df_i['sin_hs']                                                                         #水平面 H
    df_i['Id_D'] = df_i['Id'] * 0.5                                                                                                         #拡散   D
    df_i.loc[df_i['hs'] > 0, 'Id_R'] = (df_i['Id'] + df_i['Ib']) * df_i['sin_hs'] * 0.5 * 0.1                                               #反射   R

    eta = lambda c: + 2.3920 * c -3.8636 * c * c * c + 3.7568 * c * c * c * c * c - 1.3968 * c * c * c * c * c * c * c 

    df_i = df_i.fillna(0)

    df_i['Ib_E_g'] = df_i['Ib_E'] * eta(-1 * df_i['cos_hs'] * df_i['sin_AZs'])                                                              #東面   E
    df_i['Ib_S_g'] = df_i['Ib_S'] * eta(     df_i['cos_hs'] * df_i['cos_AZs'])                                                              #南面   S
    df_i['Ib_W_g'] = df_i['Ib_W'] * eta(     df_i['cos_hs'] * df_i['sin_AZs'])                                                              #西面   W
    df_i['Ib_N_g'] = df_i['Ib_N'] * eta(-1 * df_i['cos_hs'] * df_i['cos_AZs'])                                                              #北面   N
    df_i['Ib_H_g'] = df_i['Ib_H'] * eta(df_i['sin_hs'])                                                                                     #水平面 H
    df_i['Id_D_g'] = df_i['Id_D'] * 0.808                                                                                                   #拡散   D
    df_i['Id_R_g'] = df_i['Id_R'] * 0.808                                                                                                   #反射   R

    df_i = df_i.fillna(0)

    df_i['Ins_W_E'] = df_i['Ib_E']   + df_i['Id_D']   + df_i['Id_R']
    df_i['Ins_G_E'] = df_i['Ib_E_g'] + df_i['Id_D_g'] + df_i['Id_R_g']
    df_i['Ins_W_S'] = df_i['Ib_S']   + df_i['Id_D']   + df_i['Id_R']
    df_i['Ins_G_S'] = df_i['Ib_S_g'] + df_i['Id_D_g'] + df_i['Id_R_g']
    df_i['Ins_W_W'] = df_i['Ib_W']   + df_i['Id_D']   + df_i['Id_R']
    df_i['Ins_G_W'] = df_i['Ib_W_g'] + df_i['Id_D_g'] + df_i['Id_R_g']
    df_i['Ins_W_N'] = df_i['Ib_N']   + df_i['Id_D']   + df_i['Id_R']
    df_i['Ins_G_N'] = df_i['Ib_N_g'] + df_i['Id_D_g'] + df_i['Id_R_g']
    df_i['Ins_W_H'] = df_i['Ib_H']   + df_i['Id_D']   + df_i['Id_R']
    df_i['Ins_G_H'] = df_i['Ib_H_g'] + df_i['Id_D_g'] + df_i['Id_R_g']

    return df_i.fillna(0)

def make_solar(**kwargs):
    lat = kwargs['lat'] if 'lat' in kwargs else 36.00
    lon = kwargs['lon'] if 'lon' in kwargs else 140.00 
    #td  = kwargs['td']  if 'td'  in kwargs else -0.5

    if 's_ig' in kwargs:    
        s_ig = kwargs['s_ig']
        if 'td' in kwargs:
            td = kwargs['td']
        else:
            td = (s_ig.index[1] - s_ig.index[0]).seconds + (s_ig.index[1] - s_ig.index[0]).microseconds / 1000000   #t_stepの読み込み
            td = - td / 2 / 3600
        df_i = pd.concat([s_ig, sun_loc(s_ig.index, lat = lat, lon = lon, td = td)], axis = 1)
        df_i = pd.concat([df_i, sep_direct_diffuse(s_ig, df_i['hs'])], axis = 1)                                #直散分離結果の追加  
        df_i = direc_solar(df_i['Ib'], df_i['Id'],                                                              #方位別日射量の追加
                       df_i['sin_hs'], df_i['cos_hs'], df_i['hs'], df_i['sin_AZs'], df_i['cos_AZs'], df_i['AZs'])
    else:
        if 's_ib' in kwargs:    s_ib = kwargs['s_ib']
        else:                   raise Exception('ERROR: 水平面拡散日射量 s_ig がありません。法線面直達日射量 s_ib もありません。')
        if 's_id' in kwargs:    s_id = kwargs['s_id']
        else:                   raise Exception('ERROR: 水平面拡散日射量 s_ig がありません。全天日射量 s_id もありません。')
        if 'td' in kwargs:
            td = kwargs['td']
        else:
            td = (s_ib.index[1] - s_ib.index[0]).seconds + (s_ib.index[1] - s_ib.index[0]).microseconds / 1000000   #t_stepの読み込み
            td = - td / 2 / 3600
        df_i = pd.concat([s_ib, s_id, sun_loc(s_ib.index, lat = lat, lon = lon, td = td)], axis = 1)  
        df_i = direc_solar(s_ib, s_id,                                                                              #方位別日射量の追加
                       df_i['sin_hs'], df_i['cos_hs'], df_i['hs'], df_i['sin_AZs'], df_i['cos_AZs'], df_i['AZs'])

    solar =     {
        '日射_壁_E':        df_i['Ins_W_E'],
        '日射_ガラス_E':    df_i['Ins_G_E'],
        '日射_壁_S':        df_i['Ins_W_S'],
        '日射_ガラス_S':    df_i['Ins_G_S'],
        '日射_壁_W':        df_i['Ins_W_W'],
        '日射_ガラス_W':    df_i['Ins_G_W'],
        '日射_壁_N':        df_i['Ins_W_N'],
        '日射_ガラス_N':    df_i['Ins_G_N'],
        '日射_壁_H':        df_i['Ins_W_H'],
        '日射_ガラス_H':    df_i['Ins_G_H']
    }

    return df_i, solar

############################################################################################################################
#PMV PPD
############################################################################################################################
calc_R  = lambda f_cl, t_cl, t_r:               3.96e-8 * f_cl * (math.pow(t_cl + 273, 4) - math.pow(t_r + 273, 4))
calc_C  = lambda f_cl, h_c, t_cl, t_a:          f_cl * h_c * (t_cl - t_a)
calc_RC = lambda f_cl, h_c, t_cl, t_a, t_r:     calc_R(f_cl, t_cl, t_r) + calc_C(f_cl, h_c, t_cl, t_a)

def calc_PMV(Met = 1.0, W = 0.0, Clo = 1.0, t_a = 20, h_a = 50, t_r = 20, v_a = 0.2):
    M, I_cl      = Met * 58.2, Clo * 0.155
    f_cl         = (1.00 + 1.290 * I_cl) if I_cl < 0.078 else (1.05 + 0.645 * I_cl)
    t_cl         = t_a
    omega, error = 0.5, 1e12

    while abs(error) > 1e-6:
        h_c      = max(2.38 * math.pow(abs(t_cl - t_a), 0.25), 12.1 * math.sqrt(v_a))
        New_t_cl = 35.7 - 0.028 * (M - W) - I_cl * calc_RC(f_cl, h_c, t_cl, t_a, t_r)
        error    = New_t_cl - t_cl
        t_cl     = t_cl + error * omega

    E_d  = 3.05e-3 * (5733 - 6.99 * (M - W) - e(t_a, h_a))
    E_s  = 0.42 * ((M - W) - 58.15)
    E_re = 1.7e-5 * M * (5867 - e(t_a, h_a))
    C_re = 0.0014 * M * (34 - t_a)
    L    = (M - W) - E_d - E_s - E_re - C_re - calc_RC(f_cl, h_c, t_cl, t_a, t_r)
    PMV  = (0.303 * math.exp(-0.036 * M) + 0.028) * L

    return PMV

def calc_PPD(Met = 1.0, W = 0.0, Clo = 1.0, t_a = 20, h_a = 50, t_r = 20, v_a = 0.2):
    PMV = calc_PMV(Met, W, Clo, t_a, h_a, t_r, v_a)
    PPD = 100 - 95 * math.exp(-0.03353 * math.pow(PMV, 4) - 0.2179 * math.pow(PMV, 2))
    return PPD

############################################################################################################################
#Fungal Index
############################################################################################################################
def calc_fungal_index(h, t):
    a  = - 0.3
    b  = 0.685
    c1 = 0.95
    c2 = 0.07
    c3 = 25
    c4 = 7.2 
    
    x = (h - c1) / c2
    y = (t - c3) / c4

    FI = 187.25 * np.exp((((x ** 2) - 2 * a * x * y + (y ** 2)) ** b) /(2 * (a ** 2) - 2)) - 8.25

    return FI
