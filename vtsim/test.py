import vtsim as vt
import archenvlib as lib

input = {
    'index': vt.index('1h', 3600 * 10),
    'sn': {
        'Supply':   {'t_flag': vt.SN_CALC, 'x_flag': vt.SN_CALC},
        'Room1':    {'t_flag': vt.SN_CALC, 'x_flag': vt.SN_CALC, 'v': 10.0},
        'Room2':    {'t_flag': vt.SN_CALC, 'x_flag': vt.SN_CALC, 'v': 5.0},
        'Return':   {'t_flag': vt.SN_CALC, 'x_flag': vt.SN_CALC},
        'Outside':  {'t_flag': vt.SN_FIX,  'x_flag': vt.SN_FIX,  't':  5.0, 'x': lib.x(5.0, 60.0)},
    },
    'vn': {
        'Supply -> Room1 -> Return':   {'vol':  500 / 3600},
        'Supply -> Room2 -> Return':   {'vol':  500 / 3600},
        'Return -> Outside -> Supply': {'vol':  100 / 3600},
    },
    'tn': {
        'Room1 -> Outside':       {'cdtc': 100 * 0.87},
        'Room2 -> Outside':       {'cdtc': 100 * 0.87},
    },
    'aircon': {

        'AC-RAC':  {'out': 'Supply', 'set': 'Return', 'ac_model': vt.AC_RAC,
                    'ac_mode': vt.AC_HEATING, 'vol': 1000 / 3600, 
                    'pre_tmp': 22.0, 'pre_x': lib.x(22.0, 60.0),
                    'q_rtd_C': 5600.0, 'q_max_C':  5944.619999999999, 'e_rtd_C': 3.2431999999999994,
                    'q_rtd_H': 6685.3, 'q_max_H': 10047.047813999998, 'e_rtd_H': 4.157264,
                    'dualcompressor': False, 'To': 5.0, 'ho': 60.0},

#        'AC_DUCT_C': {'out': 'Supply', 'set': 'Return', 'ac_model': vt.AC_DUCT_C,
#                      'ac_mode': vt.AC_COOLING, 'vol': 1000 / 3600,
#                      'pre_tmp':       28.0,    'pre_x': lib.x(28.0, 50.0),
#                      'q_hs_rtd_H':  7733.9925600000015,  'q_hs_mid_H':  3866.9962800000008,    'q_hs_min_H':  2706.8973960000003, 
#                      'P_hs_rtd_H':  2056.9129148936177,  'V_fan_rtd_H': 1654.226845584 / 3600, 'V_fan_mid_H': 1262.113422792 / 3600, 
#                      'P_fan_rtd_H':  241.2635794112,     'P_fan_mid_H':  188.9817897056,       'V_hs_dsgn_H': 1306.8392080113601 / 3600,
#                      'q_hs_rtd_C':  7664.646360000001,   'q_hs_mid_C':  3832.3231800000003,    'q_hs_min_C':  2682.626226, 
#                      'P_hs_rtd_C':  2417.8695141955836,  'V_fan_rtd_C': 1647.195140904 / 3600, 'V_fan_mid_C': 1258.597570452 / 3600, 
#                      'P_fan_rtd_C':  240.32601878719998, 'P_fan_mid_C': 188.5130093936,        'V_hs_dsgn_C': 1301.28416131416 / 3600, 
#                      'To': 35.0,                         'ho': 60.0}

    },
    'humidity_source':{
        'h_source01':         {'set': 'Room1', 'mx': 0.0001}
    }
}

vt.run_calc(input)