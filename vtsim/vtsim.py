###############################################################################
# import
###############################################################################
import numpy as np
import pandas as pd
import time
import json
import itertools

from datetime import datetime, timedelta
import matplotlib.pyplot as plt

import vtsimc as vt

###############################################################################
# logger
###############################################################################
import logging
from pytz import timezone
from datetime import datetime

# loggerに命名する. この名前で呼び出すことで他のモジュールにも以下の設定が引き継がれる.
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

'''
# コンソールに出力するハンドラの設定
sh = logging.StreamHandler()
sh.setLevel(logging.INFO)

logger.addHandler(sh)

def customTime(*args):
    return datetime.now(timezone('Asia/Tokyo')).timetuple()

formatter = logging.Formatter(
    #fmt='%(levelname)s : %(asctime)s : %(message)s',
    fmt='%(message)s',
    datefmt="%Y-%m-%d %H:%M:%S"
)

formatter.converter = customTime
sh.setFormatter(formatter)
'''

###############################################################################
# define const
###############################################################################
SOLVE_LU:   int = vt.SOLVE_LU
SOLVE_SOR:  int = vt.SOLVE_SOR

OPT_DF:     int  = 0              #DataFrameを出力
OPT_CSV:    int  = 1              #上記に加え、csvファイルを出力
OPT_GRAPH:  int  = 2              #上記に加えグラフを描画

SN_NONE:    int = vt.SN_NONE
SN_CALC:    int = vt.SN_CALC
SN_FIX:     int = vt.SN_FIX
SN_DLY:     int = vt.SN_DLY
    
VN_SIMPLE:  int = vt.VN_SIMPLE
VN_GAP:     int = vt.VN_GAP
VN_FIX:     int = vt.VN_FIX
VN_AIRCON:  int = vt.VN_AIRCON
VN_FAN:     int = vt.VN_FAN
    
TN_SIMPLE:  int = vt.TN_SIMPLE
TN_AIRCON:  int = vt.TN_AIRCON
TN_SOLAR:   int = vt.TN_SOLAR
TN_GROUND:  int = vt.TN_GROUND
TN_HEATER:  int = vt.TN_HEATER

AC_NONE:    int = vt.AC_NONE
AC_RAC:     int = vt.AC_RAC
AC_DUCT_C:  int = vt.AC_DUCT_C
AC_RAC2:    int = vt.AC_RAC2
AC_CRIEPI:  int = vt.AC_CRIEPI

AC_AUTO:    int = vt.AC_AUTO
AC_HEATING: int = vt.AC_HEATING
AC_COOLING: int = vt.AC_COOLING
AC_STOP:    int = vt.AC_STOP

calc = vt.VTSim()

df_p   = pd.DataFrame()
df_c   = pd.DataFrame()
df_t   = pd.DataFrame()
df_x   = pd.DataFrame()
df_qv  = pd.DataFrame()
df_qt1 = pd.DataFrame()
df_qt2 = pd.DataFrame()

###############################################################################
# define lambda etc.
###############################################################################
read_csv = lambda fn:   pd.read_csv(fn, index_col = 0, parse_dates = True).fillna(method = 'bfill')\
                                                                          .fillna(method = 'ffill')     #csvファイルの読み込み

index    = lambda freq, length:   pd.DataFrame(index = pd.date_range(datetime(2021, 1, 1, 0, 0, 0), 
                                               datetime(2021, 1, 1, 0, 0, 0) + timedelta(seconds = length - 1), 
                                               freq = freq)).index                                      #頻度freq、長さlengthのindex

###############################################################################
# define function
###############################################################################
def encode(object):                                                                                     #JSON形式に変換する際のエンコード
    if isinstance(object, pd.core.indexes.datetimes.DatetimeIndex): 
        return(object.strftime('%Y/%m/%d %H:%M:%S').to_list())
    if isinstance(object, pd.core.series.Series):                   
        return(object.to_list())

def read_json(fn):                                                                                      #JSON形式を読み込みdict型で返す
    with open(fn, encoding = "utf-8") as f:
        input = json.load(f)
    return input

def read_hasp(fn):
    df = pd.DataFrame()
    str_dat = [''] * 7
    clm = ['t_ex', 'h_ex', 'i_b', 'i_d', 'n_r', 'w_d', 'w_s']

    with open(fn, 'rb') as f:
        dat = [s.decode() for s in f.readlines()]

    for d in range(365):    
        for i in range(7): 
            str_dat[i] += dat[d * 7 + i][:72]

    for c in range(7):
        df[clm[c]] = [int(str_dat[c][i * 3 : i * 3 + 3]) for i in range(24 * 365)]

    df['t_ex'] = (df['t_ex'] - 500) / 10  #気温             Exteria Temperature       ℃
    df['h_ex'] =  df['h_ex'] / 10         #重量絶対湿度     Exteria Humidity          g/kg'
    df['i_b']  =  df['i_b'] * 1.16222     #法線面直達日射   Direct Solar Insolation   kcal/(m2・h) -> W/m2 
    df['i_d']  =  df['i_d'] * 1.16222     #水平面拡散日射量 Diffuse Solar Insolation  kcal/(m2・h) -> W/m2
    df['n_r']  =  df['n_r'] * 1.16222     #夜間放射         Nocturnal Radiation       kcal/(m2・h) -> W/m2
    df['w_d'].astype(int)                 #風向             Wind Direction
                                          #0:無風,  1:NNE,  2:NE,   3:ENE, 4:E,    5:ESE,  6:SE,   7:SSE,  8:S
                                          #9:SSW,  10:SW,  11:WSW, 12:W,  13:WNW, 14:NW,  15:NNW, 16:N
    df['w_s']  = df['w_s'] / 10           #風速             Wind Speed                m/s

    return df

def write_json(input, fn):                                                                              #dict型をJSONファイルに書き出し
    input['version'] = '4.0.0'
    with open(fn, 'w') as f:
        json.dump(input, f, default = encode, ensure_ascii = False, indent = 4)

def to_json(input):                                                                                     #dict型をJSON形式に変換
    input['version'] = '4.0.0'
    return(json.dumps(input, default = encode, ensure_ascii = False, indent = 4))                       

def to_list_f(v):
    if   type(v) == list:                   return(v)
    elif type(v) == np.ndarray:             return(v)
    elif type(v) == pd.core.series.Series:  return(np.array(v))
    else:                                   return[float(v)] * calc.sts.length

def to_list_i(v):
    if   type(v) == list:                   return(v)
    elif type(v) == np.ndarray:             return(v)
    elif type(v) == pd.core.series.Series:  return(np.array(v))
    else:                                   return[int(v)] * calc.sts.length

###############################################################################
# main
###############################################################################
def run_calc(input):                                                                        #はじめに呼び出される関数    
    if   type(input) == dict:   input = to_json(input)                                      #辞書型であれば、JSON形式に変換
    elif type(input) != str:    raise Exception('ERROR: inputは、辞書型かJSON形式である必要がります。') #文字列（JSON形式)で無ければエラー
    
    input = json.loads(input)                                                               #JSON形式を辞書型に変換
    if 'vn' not in input:       input['vn'] = {}
    if 'tn' not in input:       input['tn'] = {}

    
    if 'index' in input:        set_calc_status(input)                                      #計算条件を設定
    else:                       raise Exception('ERROR: index が存在しません。')             #indexが無ければエラー

    if 'wind_pressure' in input:    input = set_wind_pressure(input)                        #外部風圧の設定
    if 'outside_temp'  in input:    input = set_outside_temp(input)                         #外気温の設定
    if 'solar'         in input:    input = set_solar(input)                                #日射量の設定
    if 'nocturnal'     in input:    input = set_nocturnal(input)                            #夜間放射の設定

    if 'room'          in input:    input = set_room(input)                                 #部屋の設定
    if 'ground'        in input:    input = set_ground(input)                               #地面の設定
    if 'wall'          in input:    input = set_wall(input)                                 #壁面の設定
    if 'glass'         in input:    input = set_glass(input)                                #ガラスの設定
    if 'room'          in input:    input = set_radiation(input)                            #放射熱伝達回路の設定

    if 'sn' in input:           input = add_capa(input)                                     #熱容量を設定
    else:                       raise Exception('ERROR: ノード(sn)が存在しません。')          #sn（ノード）が無ければエラー

    if 'aircon' in input:       input = set_aircon1(input)                                  #エアコンをセット（前処理）
    if 'heater' in input:       input = set_heater(input)                                   #ヒーターをセット

    if 'dust_source' in input:  input = set_dust_source(input)                              #発塵源のセット
    if 'air_cleaner' in input:  input = set_air_cleaner(input)                              #空気清浄機のセット

    if 'humidity_source' in input:  input = set_humidity_source(input)                      #発湿源のセット

    #with open('calc.json', 'w') as f:                                                      #計算入力を　calc.jsonに格納
    #    json.dump(input, f, ensure_ascii = False, indent = 4)

    set_sim_node(input['sn'])                                                               #sn（ノード）の設定
    set_vent_net(input['vn'])                                                               #vn（換気回路網）の設定
    set_thrm_net(input['tn'])                                                               #tn（熱回路網）の設定
    
    if 'aircon' in input:   set_aircon2(input)                                              #エアコンﾝをセット（後処理）

    logger.info('******* check vtsim c++ data *******')
    logger.info('sts     ' + str([calc.sts.length, calc.sts.t_step, calc.sts.solve, 
                                  calc.sts.step_p, calc.sts.vent_err, calc.sts.step_t, calc.sts.thrm_err,
                                  calc.sts.conv_err, calc.sts.sor_ratio, calc.sts.sor_err]))
    logger.info('sn      ' + str([n.i for n in calc.sn]))
    logger.info('node    ' + str(calc.node))
    logger.info('v_net   ' + str(calc.v_net))
    logger.info('t_net   ' + str(calc.t_net))
    logger.info('v_idc   ' + str(calc.v_idc))
    logger.info('c_idc   ' + str(calc.c_idc))
    logger.info('t_idc   ' + str(calc.t_idc))
    logger.info('x_idc   ' + str(calc.x_idc))
    logger.info('i_vn_ac ' + str(calc.i_vn_ac))
    logger.info('i_tn_ac ' + str(calc.i_tn_ac))

    logger.info('******* start vtsim c++ calc *******')
    s_time = time.time()
    calc.calc()                                                                     #計算
    e_time = time.time() - s_time    
    logger.info("calc time = {0}".format(e_time * 1000) + "[ms]")

    logger.info('******* output vtsim c++ calc *******')
    dat_list = make_df(calc.result(), pd.to_datetime(input['index'], format='%Y/%m/%d %H:%M:%S'))
    opt = input['opt'] if 'opt' in input else OPT_GRAPH

    output_calc(dat_list, opt)

def set_calc_status(input):
    logger.info('Set calc status.')
    sts  = vt.CalcStatus()

    ix = pd.to_datetime(input['index'])
    sts.length = len(ix)
    sts.t_step = (ix[1] - ix[0]).seconds + (ix[1] - ix[0]).microseconds / 1000000   #t_stepの読み込み    

    if 'calc_sts' in input: calc_sts = input['calc_sts']
    else:                   calc_sts = {}
    
    if 'solve'     in calc_sts:   sts.solve     = calc_sts['solve']
    if 'step_p'    in calc_sts:   sts.step_p    = calc_sts['step_p']
    if 'vent_err'  in calc_sts:   sts.vent_err  = calc_sts['vent_err']
    if 'step_t'    in calc_sts:   sts.step_t    = calc_sts['step_t']
    if 'thrm_err'  in calc_sts:   sts.thrm_err  = calc_sts['thrm_err']
    if 'conv_err'  in calc_sts:   sts.conv_err  = calc_sts['conv_err']
    if 'sor_ratio' in calc_sts:   sts.sor_ratio = calc_sts['sor_ratio']
    if 'sor_err'   in calc_sts:   sts.sor_err   = calc_sts['sor_err']

    calc.setup(sts)

def set_wind_pressure(input):
    logger.info('Set Wind Pressure.')
    name = ['風圧_E', '風圧_S', '風圧_W', '風圧_N', '風圧_H']

    for w in input['wind_pressure']:
        wp = input['wind_pressure']
        for nn in name:
            if nn in w: 
                nnn = nn.replace('風圧', '外部')
                if nnn not in input['sn']:  input['sn'][nnn] = {}
                input['sn'][nnn]['v_flag'] = SN_FIX
                input['sn'][nnn]['p']      =  wp

    return input

def set_outside_temp(input):
    logger.info('Set Outside Temp.')
    name = ['外部', '外部_E', '外部_S', '外部_W', '外部_N', '外部_H', '地盤', '地盤2', '地盤3']

    for t in input['outside_temp']:
        tp = input['outside_temp'][t]
        for nn in name:
            if nn == t: 
                if nn not in input['sn']:   input['sn'][nn] = {}
                input['sn'][nn] = {'t_flag': SN_FIX, 't': tp}

    return input

def set_solar(input):
    logger.info('Set Solar.')
    name = ['日射_全天_H', '日射_壁_E', '日射_壁_S', '日射_壁_W', '日射_壁_N', '日射_壁_H',
                           '日射_ガラス_E', '日射_ガラス_S', '日射_ガラス_W', '日射_ガラス_N', '日射_ガラス_H'] 
    
    for s in input['solar']:
        sl = input['solar'][s]
        for nn in name:
            if nn in s: input['sn'][nn] = {'insolation': sl}

    return input

def set_nocturnal(input):
    logger.info('Set Nocturnal Radiation.')
    name = ['夜間_H', '夜間_V']
    for n in input['nocturnal']:
        nr = input['nocturnal'][n]
        for nn in name:
            if nn in nr: input['sn'][nn] = {'insolation': nr}
    return input

def set_room(input):
    logger.info('Set Room.')
    for r in input['room']:
        rm = input['room'][r]
        if 'ground' in rm:
            if 'ground' not in input:   input['ground'] = {}
            for g in rm['ground']:
                gnd = rm['ground'][g]
                input['ground'][r + ' -> ' + g] = gnd
        if 'wall' in rm:
            if 'wall' not in input:     input['wall'] = {}
            for w in rm['wall']:
                wl = rm['wall'][w]
                input['wall'][r + '->' + w] = wl
        if 'glass' in rm:
            if 'glass' not in input:     input['glass'] = {}
            for g in rm['glass']:
                gl = rm['glass'][g]
                input['glass'][r + '->' + g] = gl
    return input

def set_ground(input):
    logger.info('Set Ground.')
    for g in input['ground']:
        
        gnd = input['ground'][g]
        n1, n2, n1_is, _ = get_n1n2(g)
        area = gnd['area']

        input['sn'][n1_is] = {'t_flag': vt.SN_CALC, 'area': area}
        input['tn'][n1 +    ' -> ' + n1_is] = {'cdtc':  area / gnd['rg']}
        input['tn'][n1_is + ' -> ' + n2   ] = {'area':  area,         'phi_0':   gnd['phi_0'],
                                               'cof_r': gnd['cof_r'], 'cof_phi': gnd['cof_phi']}
    return input

def set_wall(input):
    logger.info('Set Wall')

    for w in input['wall']:
        wl = input['wall'][w]

        n1, n2, n1_is, n1_os = get_n1n2(w)
        n1_is = n1_is + '_w'
        n1_os = n1_os + '_w'

        area    = wl['area']
        alpha_1 = wl['alpha_1'] if 'alpha_1' in wl else 4.4
        alpha_2 = wl['alpha_2'] if 'alpha_2' in wl else 23                      #外部は、まだ夜間放射考慮していない値

        input['sn'][n1_is] = {'t_flag': vt.SN_CALC, 'area': area}
        input['sn'][n1_os] = {'t_flag': vt.SN_CALC, 'area': area}

        if 't' in input['sn'][n1]:
            input['sn'][n1_is]['t'] =  input['sn'][n1]['t']
        if 't' in input['sn'][n2]:
            input['sn'][n1_os]['t'] =  input['sn'][n2]['t']

        if 'capa_w' in wl:
            input['sn'][n1_is]['capa'] =  area * wl['capa_w'] / 2
            input['sn'][n1_os]['capa'] =  area * wl['capa_w'] / 2

        input['tn'][n1    + ' -> ' + n1_is] = {'cdtc': area * alpha_1}
        input['tn'][n1_is + ' -> ' + n1_os] = {'cdtc': area * wl['U_w']}
        input['tn'][n1_os + ' -> ' + n2   ] = {'cdtc': area * alpha_2}

        if 'solar' in wl:
            input['tn'][n1_os + ' -> ' + wl['solar']] = {'ms': area * wl['eta_w']}

        if 'nocturnal' in wl:
            input['tn'][n1_os + ' -> ' + wl['nocturnal']] = {'ms': area * wl['epsilon_w']}

    return input

def set_glass(input):
    logger.info('Set Glass')
    for g in input['glass']:
        gl = input['glass'][g]

        n1, n2, n1_is, n1_os = get_n1n2(g)
        n1_is = n1_is + '_g'
        n1_os = n1_os + '_g'

        area    = gl['area']
        alpha_1 = gl['alpha_1'] if 'alpha_1' in gl else 9.0
        alpha_2 = gl['alpha_2'] if 'alpha_2' in gl else 25.0

        input['sn'][n1_is] = {'t_flag': vt.SN_CALC, 'area': area}
        input['sn'][n1_os] = {'t_flag': vt.SN_CALC, 'area': area}

        if 't' in input['sn'][n1]:
            input['sn'][n1_is]['t'] =  input['sn'][n1]['t']
        if 't' in input['sn'][n2]:
            input['sn'][n1_os]['t'] =  input['sn'][n2]['t']

        input['tn'][n1    + ' -> ' + n1_is] = {'cdtc': area * alpha_1}
        input['tn'][n1_is + ' -> ' + n1_os] = {'cdtc': area * gl['U_w']}
        input['tn'][n1_os + ' -> ' + n2   ] = {'cdtc': area * alpha_2}

        if 'solar' in gl:
            input['tn'][n1 + ' -> ' + gl['solar']] = {'ms': area * gl['eta_w']}

        if 'nocturnal' in gl:
            input['tn'][n1_os + ' -> ' + gl['nocturnal']] = {'ms': area * gl['epsilon_w']}

    return input

def set_radiation(input):
    logger.info('Set Radiation.')
    for r in input['room']:
        node = [n for n in input['sn'] if n.find(r + '_') == 0 and '_is' in n]
        sum_area = sum([input['sn'][n]['area'] for n in node])
        for pair in itertools.combinations(node, 2):
            area0 = input['sn'][pair[0]]['area']
            area1 = input['sn'][pair[1]]['area']
            input['tn'][pair[0] + ' -> ' + pair[1]] = {'cdtc': area0 * area1 / sum_area * 4.6}

    return input

def add_capa(input):    
    logger.info('Add Capacity.')
    for n in [n for n in input['sn'] if 'capa' in input['sn'][n]]:                              #熱容量の設定のあるノード
        nc = n + '_c'
        
        input['sn'][nc] = {'t_flag': vt.SN_DLY, 's_i': n}                                       #計算フラグ、親ノードの設定
        if 't' in input['sn'][n]:   input['sn'][nc]['t'] = input['sn'][n]['t']                  #初期温度の継承

        input['tn'][n + ' -> ' + nc] = {'cdtc': input['sn'][n]['capa'] / calc.sts.t_step}       #熱容量の設定 

    return input

def set_aircon1(input):
    logger.info('Set Aircon (Pre-processing).')
    for a in input['aircon']:
        ac = input['aircon'][a]
        ac_in, ac_out = a + '_in', a + '_out'

        if 'set' in ac:             n3 = ac['set']
        else:                       raise Exception('ERROR: エアコンのsetが設定されていません')
        n1 = ac['in']  if 'in'  in ac else ac['set']
        n2 = ac['out'] if 'out' in ac else ac['set']
        
        if 'x_flag' in input['sn'][n1]:
            if  input['sn'][n1]['x_flag'] == vt.SN_CALC:
                input['sn'][ac_in]            = {'t_flag': vt.SN_CALC, 'x_flag': vt.SN_CALC}
                input['sn'][ac_out]           = {'t_flag': vt.SN_CALC, 'x_flag': vt.SN_CALC}
            else:
                input['sn'][ac_in]            = {'t_flag': vt.SN_CALC}
                input['sn'][ac_out]           = {'t_flag': vt.SN_CALC}
        else:
            input['sn'][ac_in]            = {'t_flag': vt.SN_CALC}
            input['sn'][ac_out]           = {'t_flag': vt.SN_CALC}

        vol = ac['vol'] if 'vol' in ac else 1000 / 3600

        input['vn'][n1     + ' -> ' + ac_in]  = {'vol': vol}
        input['vn'][ac_in  + ' -> ' + ac_out] = {'ac_vol': vol}
        input['vn'][ac_out + ' -> ' + n2]     = {'vol': vol}

        if ac['ac_model'] == AC_NONE:
            input['tn'][n3  + ' -> ' + ac_out] = {'ac_model': ac['ac_model'], 'ac_mode': ac['ac_mode'], 'pre_tmp': ac['pre_tmp']}

        if ac['ac_model'] == AC_RAC:
            input['tn'][n3  + ' -> ' + ac_out] = {'ac_model': ac['ac_model'], 'ac_mode': ac['ac_mode'], 
                                                  'pre_tmp': ac['pre_tmp'], 'pre_x': ac['pre_x'],
                                                  'q_rtd_C': ac['q_rtd_C'], 'q_max_C': ac['q_max_C'], 'e_rtd_C': ac['e_rtd_C'],
                                                  'q_rtd_H': ac['q_rtd_H'], 'q_max_H': ac['q_max_H'], 'e_rtd_H': ac['e_rtd_H'],
                                                  'dualcompressor': ac['dualcompressor'], 'To': ac['To'], 'ho': ac['ho']}
        if ac['ac_model'] == AC_DUCT_C:
            input['tn'][n3  + ' -> ' + ac_out] = {'ac_model':    ac['ac_model'],    'ac_mode':     ac['ac_mode'], 
                                                  'pre_tmp':     ac['pre_tmp'],     'pre_x':       ac['pre_x'],
                                                  'q_hs_rtd_H':  ac['q_hs_rtd_H'],  'q_hs_mid_H':  ac['q_hs_mid_H'], 
                                                  'q_hs_min_H':  ac['q_hs_min_H'],  'P_hs_rtd_H':  ac['P_hs_rtd_H'],  
                                                  'V_fan_rtd_H': ac['V_fan_rtd_H'], 'V_fan_mid_H': ac['V_fan_mid_H'], 
                                                  'P_fan_rtd_H': ac['P_fan_rtd_H'], 'P_fan_mid_H': ac['P_fan_mid_H'], 
                                                  'V_hs_dsgn_H': ac['V_hs_dsgn_H'],
                                                  'q_hs_rtd_C':  ac['q_hs_rtd_C'],  'q_hs_mid_C':  ac['q_hs_mid_C'], 
                                                  'q_hs_min_C':  ac['q_hs_min_C'],  'P_hs_rtd_C':  ac['P_hs_rtd_C'], 
                                                  'V_fan_rtd_C': ac['V_fan_rtd_C'], 'V_fan_mid_C': ac['V_fan_mid_C'], 
                                                  'P_fan_rtd_C': ac['P_fan_rtd_C'], 'P_fan_mid_C': ac['P_fan_mid_C'], 
                                                  'V_hs_dsgn_C': ac['V_hs_dsgn_C'], 
                                                  'To':          ac['To'],           'ho':         ac['ho']}

    return input

def set_aircon2(input):
    logger.info('Set Aircon (Post-processing).')
    for a in input['aircon']:
        ac = input['aircon'][a]
        ac_in, ac_out, n3 = a + '_in', a + '_out', ac['set']
        
        for i, nt in enumerate(calc.vn):
            if nt.vn_type == vt.VN_AIRCON:
                if (calc.node[ac_in] == nt.i1) and (calc.node[ac_out] == nt.i2):    calc.vn_aircon_add(i)
        for i, nt in enumerate(calc.tn):
            if nt.tn_type == vt.TN_AIRCON:
                if (calc.node[n3]    == nt.i1) and (calc.node[ac_out] == nt.i2):    calc.tn_aircon_add(i)

def set_heater(input):
    logger.info('Set Heater.')
    for h in input['heater']:
        ht = input['heater'][h]
        input['sn'][h] = {'h_input': ht['output']}
        input['tn'][ht['set'] + ' -> ' + h] = {}

    return input

def set_dust_source(input):
    logger.info('Set Dust Source.')
    for d in input['dust_source']:
        ds = input['dust_source'][d]
        input['sn'][ds['set']]['m'] = ds['m']
    return input

def set_air_cleaner(input):
    logger.info('Set Air Cleaner.')
    for a in input['air_cleaner']:
        ac = input['air_cleaner'][a]
        input['vn'][ac['set']]['eta'] = ac['eta']
    return input

def set_humidity_source(input):
    logger.info('Set Humidity Source.')
    for m in input['humidity_source']:
        ms = input['humidity_source'][m]
        input['sn'][ms['set']]['mx'] = ms['mx']
    return input

def sep_sfx(s, opt = False):
    if s.find(':') != -1:       
        if opt: 
            return s[:s.find(':')] + '(' + s[s.find(':') + 1:] + ')'
        else:
            print(s + ' の添字 ' + s[s.find(':'):] + ' は無視されました。')
            return s[:s.find(':')]
    else:                      
        return s

def get_node(s):
    o_node = []
    for n in s.replace(' ','').split(','):  o_node.append(sep_sfx(n))
    return o_node

def set_sim_node(sn):
    logger.info('Set SimNode.')
    i = 0
    for n in sn:
        for nn in get_node(n):
            logger.info('    setting node ' + str(i) + ' : ' + nn)
            calc.set_node(nn, i)
            v_flag = sn[n]['v_flag'] if 'v_flag' in sn[n] else vt.SN_NONE
            c_flag = sn[n]['c_flag'] if 'c_flag' in sn[n] else vt.SN_NONE
            t_flag = sn[n]['t_flag'] if 't_flag' in sn[n] else vt.SN_NONE
            x_flag = sn[n]['x_flag'] if 'x_flag' in sn[n] else vt.SN_NONE
            calc.sn_add(i, [v_flag, c_flag, t_flag, x_flag])

            if 'p'           in sn[n]:    calc.sn[i].p     = to_list_f(sn[n]['p'])                    #圧力、行列で設定可能
            if 'c'           in sn[n]:    calc.sn[i].c     = to_list_f(sn[n]['c'])                    #濃度、行列で設定可能
            if 't'           in sn[n]:    calc.sn[i].t     = to_list_f(sn[n]['t'])                    #温度、行列で設定可能
            if 'x'           in sn[n]:    calc.sn[i].x     = to_list_f(sn[n]['x'])                    #絶対湿度、行列で設定可能  
            if 'insolation'  in sn[n]:    calc.sn[i].h_sr  = to_list_f(sn[n]['insolation'])           #日射量、行列で設定可能
            if 'h_input'     in sn[n]:    calc.sn[i].h_inp = to_list_f(sn[n]['h_input'])              #発熱、行列で設定可能
            if 'v'           in sn[n]:    calc.sn[i].v     = to_list_f(sn[n]['v'])                    #気積、行列で設定可能
            if 'm'           in sn[n]:    calc.sn[i].m     = to_list_f(sn[n]['m'])                    #発生量、行列で設定可能
            if 'beta'        in sn[n]:    calc.sn[i].beta  = to_list_f(sn[n]['beta'])                 #濃度減少率、行列で設定可能
            if 'mx'          in sn[n]:    calc.sn[i].mx    = to_list_f(sn[n]['mx'])                   #発湿量、行列で設定可能
            if 's_i'         in sn[n]:    calc.sn[i].s_i   = calc.node[sn[nn]['s_i']]

            i = i + 1

def get_n1n2(nt):  
    s = nt.replace(' ', '')

    if s.find(':') == -1:   s, sfx = s, ''
    else:                   s, sfx = s[:s.find(':')], '(' + s[s.find(':') + 1:] + ')'

    if s.find('->') == -1:  raise Exception('ERROR: vnもしくはtnのキーに -> が存在しません')
    if s.count('->') >= 2:  raise Exception('ERROR: キーに -> が2回以上出現します')

    n1 = s.split('->')[0]
    n2 = s.split('->')[1]
    n1_is = n1 + '_' + n2 + '_is' + sfx
    n1_os = n1 + '_' + n2 + '_os' + sfx

    return n1, n2, n1_is, n1_os

def get_network(s):
    o_network = []
    for nt in s.replace(' ', '').split(','):
        if nt.find('->') == -1:
            raise Exception('ERROR: vnもしくはtnのキーに -> が存在しません')
        else:
            nodes = nt.split('->') 
            for i in range(len(nodes) - 1):
                o_network.append(sep_sfx(nodes[i]) + ' -> ' + sep_sfx(nodes[i + 1]))

    return o_network

def set_vent_net(vn):
    logger.info('Set VentNet.')
    i = 0
    for nt in vn:
        for n1n2 in get_network(nt):
            logger.info('    setting ventnet ' + str(i) + ' : ' + n1n2)
            if 'type' not in vn[nt]:
                if   ('alpha'  in vn[nt]) and ('area' in vn[nt]):  vn_type = vt.VN_SIMPLE
                elif ('a'      in vn[nt]) and ('n'    in vn[nt]):  vn_type = vt.VN_GAP
                elif ('qmax'   in vn[nt]) and ('pmax' in vn[nt]):  vn_type = vt.VN_FAN
                elif  'vol'    in vn[nt]:                          vn_type = vt.VN_FIX
                elif  'ac_vol' in vn[nt]:                          vn_type = vt.VN_AIRCON
                else:                                             raise Exception('ERROR: ' + nt + 'のvn_typeを認識できません')    
            else:
                vn_type = vn[nt]['type']

            h1 = to_list_f(vn[nt]['h1']) if 'h1' in vn[nt] else to_list_f(0.0)                          #高さ1、行列設定不可
            h2 = to_list_f(vn[nt]['h2']) if 'h2' in vn[nt] else to_list_f(0.0)                          #高さ2、行列設定不可
            n1, n2, _, _ = get_n1n2(n1n2)
            if n1 not in calc.node: raise Exception('ERROR: ノード(sn)の中に ' + n1 + ' がありません。')
            if n2 not in calc.node: raise Exception('ERROR: ノード(sn)の中に ' + n2 + ' がありません。')

            nn = n1n2
            j = 1
            while n1n2 in calc.v_net.keys():
                n1n2 = nn + ' : ' + str(j)
                j = j + 1

            calc.vn_add(i, n1n2, calc.node[n1], calc.node[n2], vn_type, h1, h2)
            
            if vn_type == vt.VN_FIX:   
                calc.vn[i].qv = to_list_f(vn[nt]['vol'])                                    #風量固定値、行列で設定可能
            if vn_type == vt.VN_AIRCON:
                calc.vn[i].qv = to_list_f(vn[nt]['ac_vol'])                                 #風量固定値、行列で設定可能
            if vn_type == vt.VN_SIMPLE:                                                                 
                calc.vn[i].alpha = to_list_f(vn[nt]['alpha'])
                calc.vn[i].area  = to_list_f(vn[nt]['area'])                                #単純開口、行列で設定可能
            if vn_type == vt.VN_GAP:           
                calc.vn[i].a     = to_list_f(vn[nt]['a'])
                calc.vn[i].n     = to_list_f(vn[nt]['n'])                                   #隙間、行列で設定可能
            if vn_type == vt.VN_FAN:           
                calc.vn[i].q_max = to_list_f(vn[nt]['qmax']) 
                calc.vn[i].p_max = to_list_f(vn[nt]['pmax']) 
                calc.vn[i].q1    = to_list_f(vn[nt]['q1'])
                calc.vn[i].p1    = to_list_f(vn[nt]['p1'])                                  #ファン、行列で設定可能
            
            if 'eta' in vn[nt]:
                calc.vn[i].eta = to_list_f(vn[nt]['eta'])                                   #粉じん除去率、行列で設定可能

            i = i + 1

def set_thrm_net(tn):
    logger.info('Set ThrmNet.')
    i = 0
    for nt in tn:
        for n1n2 in get_network(nt):
            logger.info('    setting thrmnet ' + str(i) + ' : ' + n1n2)
            if 'type' not in tn[nt]:
                if    'ms'       in tn[nt]:                            tn_type = vt.TN_SOLAR
                elif ('phi_0'    in tn[nt]) and ('cof_r'   in tn[nt]): tn_type = vt.TN_GROUND
                elif ('ac_mode'  in tn[nt]) and ('pre_tmp' in tn[nt]): tn_type = vt.TN_AIRCON
                elif  'cdtc'     in tn[nt]:                            tn_type = vt.TN_SIMPLE
                elif tn[nt] == {}:                                     tn_type = vt.TN_HEATER
                else:                                                  raise Exception('ERROR: ' + nt + 'のtn_typeを認識できません') 
            else:
                tn_type = tn[nt]['type']

            n1, n2, _, _ = get_n1n2(n1n2)
            if n1 not in calc.node: raise Exception('ERROR: ノード(sn)の中に ' + n1 + ' がありません。')
            if n2 not in calc.node: raise Exception('ERROR: ノード(sn)の中に ' + n2 + ' がありません。')

            nn = n1n2
            j = 1
            while n1n2 in calc.t_net.keys():
                n1n2 = nn + ' : ' + str(j)
                j = j + 1

            calc.tn_add(i, n1n2, calc.node[n1], calc.node[n2], tn_type)

            if tn_type == vt.TN_SIMPLE:     
                calc.tn[i].cdtc = to_list_f(tn[nt]['cdtc'])                                #コンダクタンス、行列設定可能
            if tn_type == vt.TN_AIRCON:
                if tn[nt]['ac_model'] == AC_NONE:
                    calc.tn[i].ac_model       = tn[nt]['ac_model']
                    calc.tn[i].ac_mode        = to_list_i(tn[nt]['ac_mode']) 
                    calc.tn[i].pre_tmp        = to_list_f(tn[nt]['pre_tmp'])
                if tn[nt]['ac_model'] == AC_RAC:
                    calc.tn[i].ac_model       = tn[nt]['ac_model']
                    calc.tn[i].ac_mode        = to_list_i(tn[nt]['ac_mode']) 
                    calc.tn[i].pre_tmp        = to_list_f(tn[nt]['pre_tmp'])
                    calc.tn[i].pre_x          = to_list_f(tn[nt]['pre_x'])
                    calc.tn[i].q_rtd_C        = tn[nt]['q_rtd_C']                               #エアコン運転モード
                    calc.tn[i].q_max_C        = tn[nt]['q_max_C']                               #エアコン運転モード
                    calc.tn[i].e_rtd_C        = tn[nt]['e_rtd_C']                               #エアコン運転モード
                    calc.tn[i].q_rtd_H        = tn[nt]['q_rtd_H']                               #エアコン運転モード
                    calc.tn[i].q_max_H        = tn[nt]['q_max_H']                               #エアコン運転モード
                    calc.tn[i].e_rtd_H        = tn[nt]['e_rtd_H']                               #エアコン運転モード
                    calc.tn[i].dualcompressor = tn[nt]['dualcompressor']                        #エアコン運転モード
                    calc.tn[i].To             = to_list_f(tn[nt]['To'])                         #エアコン運転モード
                    calc.tn[i].ho             = to_list_f(tn[nt]['ho'])                         # エアコン運転モード
                if tn[nt]['ac_model'] == AC_DUCT_C:
                    calc.tn[i].ac_model       = tn[nt]['ac_model']
                    calc.tn[i].ac_mode        = to_list_i(tn[nt]['ac_mode']) 
                    calc.tn[i].pre_tmp        = to_list_f(tn[nt]['pre_tmp'])
                    calc.tn[i].pre_x          = to_list_f(tn[nt]['pre_x'])  
                    calc.tn[i].q_hs_rtd_H     = tn[nt]['q_hs_rtd_H']
                    calc.tn[i].q_hs_mid_H     = tn[nt]['q_hs_mid_H']
                    calc.tn[i].q_hs_min_H     = tn[nt]['q_hs_min_H'] 
                    calc.tn[i].P_hs_rtd_H     = tn[nt]['P_hs_rtd_H']
                    calc.tn[i].V_fan_rtd_H    = tn[nt]['V_fan_rtd_H'] 
                    calc.tn[i].V_fan_mid_H    = tn[nt]['V_fan_mid_H'] 
                    calc.tn[i].P_fan_rtd_H    = tn[nt]['P_fan_rtd_H'] 
                    calc.tn[i].P_fan_mid_H    = tn[nt]['P_fan_mid_H'] 
                    calc.tn[i].V_hs_dsgn_H    = tn[nt]['V_hs_dsgn_H']
                    calc.tn[i].q_hs_rtd_C     = tn[nt]['q_hs_rtd_C']  
                    calc.tn[i].q_hs_mid_C     = tn[nt]['q_hs_mid_C'] 
                    calc.tn[i].q_hs_min_C     = tn[nt]['q_hs_min_C']  
                    calc.tn[i].P_hs_rtd_C     = tn[nt]['P_hs_rtd_C']
                    calc.tn[i].V_fan_rtd_C    = tn[nt]['V_fan_rtd_C']
                    calc.tn[i].V_fan_mid_C    = tn[nt]['V_fan_mid_C']
                    calc.tn[i].P_fan_rtd_C    = tn[nt]['P_fan_rtd_C']
                    calc.tn[i].P_fan_mid_C    = tn[nt]['P_fan_mid_C']
                    calc.tn[i].V_hs_dsgn_C    = tn[nt]['V_hs_dsgn_C']
                    calc.tn[i].To             = to_list_f(tn[nt]['To'])                    #エアコン運転モード
                    calc.tn[i].ho             = to_list_f(tn[nt]['ho'])                    # エアコン運転モード

            if tn_type == vt.TN_SOLAR:       
                calc.tn[i].ms      = to_list_f(tn[nt]['ms'])                               #日射熱取得率、行列設定可能
            if tn_type == vt.TN_GROUND:     
                calc.tn[i].area    = to_list_f(tn[nt]['area'])           
                calc.tn[i].phi_0   = tn[nt]['phi_0']
                calc.tn[i].cof_r   = tn[nt]['cof_r']
                calc.tn[i].cof_phi = tn[nt]['cof_phi']                                     #地盤熱応答、行列設定不可（面積と断熱性能はOK）         
            i = i + 1

def get_keys(dict):
    keys = []
    for i in range(len(dict)):
        keys.append([k for k, v in dict.items() if v == i][0])
    return keys

def make_df(res, ix):
    global df_p, df_c, df_t, df_x, df_qv, df_qt1, df_qt2, df_elec
    dat_list = []

    if len(res[0]) != 0:    
        df_p  = pd.DataFrame(np.array(res[0]).T,  index = ix, columns = get_keys(calc.node))
        dat_list.append({'fn': 'vent_p.csv',   'title': '圧力',  'unit': '[Pa]', 'df': df_p})
    else:
        df_p = None

    if len(res[1]) != 0:    
        df_c  = pd.DataFrame(np.array(res[1]).T,  index = ix, columns = get_keys(calc.node))
        dat_list.append({'fn': 'vent_c.csv',   'title': '濃度',  'unit': '[個/L]', 'df': df_c})
    else:
        df_c = None

    if len(res[2]) != 0:    
        df_t  = pd.DataFrame(np.array(res[2]).T,  index = ix, columns = get_keys(calc.node))
        dat_list.append({'fn': 'thrm_t.csv',   'title': '温度',  'unit': '[℃]', 'df': df_t})
    else:
        df_t = None

    if len(res[3]) != 0:    
        df_x = pd.DataFrame(np.array(res[3]).T,  index = ix, columns = get_keys(calc.node))
        dat_list.append({'fn': 'vent_x.csv',   'title': '絶対湿度',  'unit': '[kg/kg(DA)]', 'df': df_x})
    else:
        df_x = None
    
    if len(res[4]) != 0:    
        df_qv  = pd.DataFrame(np.array(res[4]).T,  index = ix, columns = get_keys(calc.v_net))
        dat_list.append({'fn': 'vent_qv.csv',  'title': '風量',  'unit': '[m3/s]', 'df': df_qv})
    else:
        df_qv = None

    if len(res[5]) != 0:    
        df_qt1 = pd.DataFrame(np.array(res[5]).T,  index = ix, columns = get_keys(calc.v_net))
        dat_list.append({'fn': 'thrm_qt1.csv', 'title': '熱量1', 'unit': '[W]', 'df': df_qt1})
    else:
        df_qt1 = None

    if len(res[6]) != 0:    
        df_qt2 = pd.DataFrame(np.array(res[6]).T,  index = ix, columns = get_keys(calc.t_net))
        dat_list.append({'fn': 'thrm_qt2.csv', 'title': '熱量2', 'unit': '[W]', 'df': df_qt2})
    else:
        df_qt2 = None

    if len(res[7]) != 0:    
        df_elec = pd.DataFrame(np.array(res[7]).T,  index = ix, columns = get_keys(calc.t_net))
        dat_list.append({'fn': 'thrm_elec.csv', 'title': '電力', 'unit': '[kW/h]', 'df': df_elec})
    else:
        df_elec = None

    return dat_list

def output_calc(dat_list, opt):
    if opt == OPT_CSV or opt == OPT_GRAPH:
        logger.info('Output csv files.')
        for d in dat_list:      d['df'].to_csv(d['fn'], encoding = 'utf_8_sig')

    if opt == OPT_GRAPH:
        logger.info('Draw Graphs.')
        fig = plt.figure(facecolor = 'w', figsize = (18, len(dat_list) * 4))
        fig.subplots_adjust(wspace = -0.1, hspace=0.9)

        for i, d in enumerate(dat_list):
            a = fig.add_subplot(len(dat_list), 1, i + 1)
            for cl in d['df'].columns:
                a.plot(d['df'][cl], linewidth = 1.0, label = cl)
            a.legend(ncol = 5, bbox_to_anchor = (0, 1.05, 1, 0), 
                     loc = 'lower right', borderaxespad = 0, facecolor = 'w', edgecolor = 'k')
            a.set_title(d['title'], loc='left')
            a.set_ylabel(d['unit'])