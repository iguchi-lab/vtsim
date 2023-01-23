#暖冷房期間
period_1      = [1] * 158 + [0] * 32 + [-1] *  53 + [0] * 23 + [1] * 99
period_2      = [1] * 155 + [0] * 40 + [-1] *  48 + [0] * 25 + [1] * 97
period_3      = [1] * 151 + [0] * 39 + [-1] *  53 + [0] * 29 + [1] * 93
period_4      = [1] * 150 + [0] * 40 + [-1] *  53 + [0] * 30 + [1] * 92
period_5      = [1] * 135 + [0] * 51 + [-1] *  57 + [0] * 39 + [1] * 83
period_6      = [1] * 111 + [0] * 38 + [-1] * 117 + [0] * 41 + [1] * 58
period_7      = [1] *  86 + [0] * 48 + [-1] * 152 + [0] * 43 + [1] * 36
period_8      = [1] *   0 + [0] * 83 + [-1] * 265 + [0] * 17 + [1] * 0

#休日
holiday      = [1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,
                1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,
                1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,  1, 0, 1, 0, 0, 0, 1,
                1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,
                1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 1, 1, 1, 1,  1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,
                1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,
                1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,
                1, 0, 0, 0, 1, 0, 1,  1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,
                1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,
                1, 0, 0, 0, 0, 1, 1,  1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,
                1, 0, 1, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 1, 1,
                1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 1, 0, 1,  1, 0, 0, 0, 0, 0, 1,
                1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1,
                1]

#主たる居室の暖冷房スケジュール
ac_mode_w_MR_h  = [3, 3, 3,  3, 3, 3,  1, 1, 1,  1, 3, 3,  1, 1, 3,  3, 1, 1,  1, 1, 1,  1, 1, 1]
ac_mode_h_MR_h  = [3, 3, 3,  3, 3, 3,  3, 3, 1,  1, 1, 1,  1, 1, 3,  3, 1, 1,  1, 1, 1,  1, 1, 3]
ac_mode_w_MR_c  = [3, 3, 3,  3, 3, 3,  2, 2, 2,  2, 3, 3,  2, 2, 3,  3, 2, 2,  2, 2, 2,  2, 2, 2]
ac_mode_h_MR_c  = [3, 3, 3,  3, 3, 3,  3, 3, 2,  2, 2, 2,  2, 2, 3,  3, 2, 2,  2, 2, 2,  2, 2, 3]

#寝室の暖冷房スケジュール
ac_mode_w_BR_h  = [3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3]
ac_mode_h_BR_h  = [3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3]
ac_mode_w_BR_c  = [2, 2, 2,  2, 2, 2,  2, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 2]
ac_mode_h_BR_c  = [2, 2, 2,  2, 2, 2,  2, 2, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 2]

#子供室1の暖冷房スケジュール
ac_mode_w_CR1_h = [3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 1,  3, 1, 1]
ac_mode_h_CR1_h = [3, 3, 3,  3, 3, 3,  3, 3, 1,  1, 1, 1,  3, 3, 3,  3, 1, 1,  1, 3, 1,  1, 1, 3]
ac_mode_w_CR1_c = [2, 2, 2,  2, 2, 2,  2, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 2,  3, 2, 2]
ac_mode_h_CR1_c = [2, 2, 2,  2, 2, 2,  2, 2, 2,  2, 2, 2,  3, 3, 3,  3, 2, 2,  2, 3, 2,  2, 2, 2]

#子供室2の暖冷房スケジュール
ac_mode_w_CR2_h = [3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  1, 3, 3,  1, 1, 3]
ac_mode_h_CR2_h = [3, 3, 3,  3, 3, 3,  3, 3, 3,  1, 1, 1,  1, 3, 3,  3, 3, 3,  3, 3, 1,  1, 1, 3]
ac_mode_w_CR2_c = [2, 2, 2,  2, 2, 2,  2, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  2, 3, 3,  2, 2, 2]
ac_mode_h_CR2_c = [2, 2, 2,  2, 2, 2,  2, 2, 3,  2, 2, 2,  2, 3, 3,  3, 3, 3,  3, 3, 2,  2, 2, 2]


#主たる居室の設定温度
pre_tmp_w_MR_h  = [20.0] * 24
pre_tmp_h_MR_h  = [20.0] * 24
pre_tmp_w_MR_c  = [27.0] * 24 
pre_tmp_h_MR_c  = [27.0] * 24

#寝室の設定温度
pre_tmp_w_BR_h  = [20.0] * 24
pre_tmp_h_BR_h  = [20.0] * 24
pre_tmp_w_BR_c  = [28.0] * 24 
pre_tmp_h_BR_c  = [28.0] * 24 

#子供室1の設定温度
pre_tmp_w_CR1_h = [20.0] * 24
pre_tmp_h_CR1_h = [20.0] * 24
pre_tmp_w_CR1_c = [28.0] *  8 + [27.0] * 15 + [28.0] * 1
pre_tmp_h_CR1_c = [28.0] *  8 + [27.0] * 15 + [28.0] * 1

#子供室2の設定温度
pre_tmp_w_CR2_h = [20] * 24
pre_tmp_h_CR2_h = [20] * 24
pre_tmp_w_CR2_c = [28.0] *  8 + [27.0] * 15 + [28.0] * 1
pre_tmp_h_CR2_c = [28.0] *  8 + [27.0] * 15 + [28.0] * 1


vnt_w_MR        = [0, 0, 0,  0, 0, 0,  75, 0,  0,      0,   0, 0,    75, 0, 0,  0,    0,     0,   150,   150,    0,   0,   0,   0]
vnt_h_MR        = [0, 0, 0,  0, 0, 0,   0, 0, 75,      0,   0, 0,    75, 0, 0,  0,    0,   150,   150,     0,    0,   0,   0,   0]
vnt_w_BT        = [0, 0, 0,  0, 0, 0,   0, 0,  0,      0,   0, 0,     0, 0, 0,  0,    0,     0,     0,     0,   50,  25, 100,   0]
vnt_h_BT        = [0, 0, 0,  0, 0, 0,   0, 0,  0,      0,   0, 0,     0, 0, 0,  0,   75,    25,     0,     0,   25,  25, 100,   0]
vnt_w_RR        = [0, 0, 0,  0, 0, 6,   2, 0,  0.8,    0,   0, 0.8,   0, 0, 0,  0.8,  0.8,   0.8,   0.8,   0.8,  2,   0,   2.8, 0]
vnt_h_RR        = [0, 0, 0,  0, 0, 0,   4, 4,  0,    1.2, 1.2, 0,     0, 0, 0,  2,    0.8,   0,    -2,   0.8,    0,   2,   0.8, 0]

def make_8760_data(ac_mode, holiday, w_h, h_h, w_c, h_c, default):
    output = []
    for i in range(365):
        if  ((ac_mode[i] ==  1) and (holiday[i] == 1)):   output += h_h
        elif((ac_mode[i] ==  1) and (holiday[i] == 0)):   output += w_h
        elif (ac_mode[i] ==  0):                          output += [default] * 24
        elif((ac_mode[i] == -1) and (holiday[i] == 1)):   output += h_c
        elif((ac_mode[i] == -1) and (holiday[i] == 0)):   output += w_c
    return output

ac_mode = {'region1': {'MR':  make_8760_data(period_1, holiday, ac_mode_w_MR_h,  ac_mode_h_MR_h,  ac_mode_w_MR_c,  ac_mode_h_MR_c,  3), 
                       'BR':  make_8760_data(period_1, holiday, ac_mode_w_BR_h,  ac_mode_h_BR_h,  ac_mode_w_BR_c,  ac_mode_h_BR_c,  3), 
                       'CR1': make_8760_data(period_1, holiday, ac_mode_w_CR1_h, ac_mode_h_CR1_h, ac_mode_w_CR1_c, ac_mode_h_CR1_c, 3), 
                       'CR2': make_8760_data(period_1, holiday, ac_mode_w_CR2_h, ac_mode_h_CR2_h, ac_mode_w_CR2_c, ac_mode_h_CR2_c, 3)},
           'region2': {'MR':  make_8760_data(period_2, holiday, ac_mode_w_MR_h,  ac_mode_h_MR_h,  ac_mode_w_MR_c,  ac_mode_h_MR_c,  3), 
                       'BR':  make_8760_data(period_2, holiday, ac_mode_w_BR_h,  ac_mode_h_BR_h,  ac_mode_w_BR_c,  ac_mode_h_BR_c,  3), 
                       'CR1': make_8760_data(period_2, holiday, ac_mode_w_CR1_h, ac_mode_h_CR1_h, ac_mode_w_CR1_c, ac_mode_h_CR1_c, 3), 
                       'CR2': make_8760_data(period_2, holiday, ac_mode_w_CR2_h, ac_mode_h_CR2_h, ac_mode_w_CR2_c, ac_mode_h_CR2_c, 3)},
           'region3': {'MR':  make_8760_data(period_3, holiday, ac_mode_w_MR_h,  ac_mode_h_MR_h,  ac_mode_w_MR_c,  ac_mode_h_MR_c,  3), 
                       'BR':  make_8760_data(period_3, holiday, ac_mode_w_BR_h,  ac_mode_h_BR_h,  ac_mode_w_BR_c,  ac_mode_h_BR_c,  3), 
                       'CR1': make_8760_data(period_3, holiday, ac_mode_w_CR1_h, ac_mode_h_CR1_h, ac_mode_w_CR1_c, ac_mode_h_CR1_c, 3), 
                       'CR2': make_8760_data(period_3, holiday, ac_mode_w_CR2_h, ac_mode_h_CR2_h, ac_mode_w_CR2_c, ac_mode_h_CR2_c, 3)},
           'region4': {'MR':  make_8760_data(period_4, holiday, ac_mode_w_MR_h,  ac_mode_h_MR_h,  ac_mode_w_MR_c,  ac_mode_h_MR_c,  3), 
                       'BR':  make_8760_data(period_4, holiday, ac_mode_w_BR_h,  ac_mode_h_BR_h,  ac_mode_w_BR_c,  ac_mode_h_BR_c,  3), 
                       'CR1': make_8760_data(period_4, holiday, ac_mode_w_CR1_h, ac_mode_h_CR1_h, ac_mode_w_CR1_c, ac_mode_h_CR1_c, 3), 
                       'CR2': make_8760_data(period_4, holiday, ac_mode_w_CR2_h, ac_mode_h_CR2_h, ac_mode_w_CR2_c, ac_mode_h_CR2_c, 3)},
           'region5': {'MR':  make_8760_data(period_5, holiday, ac_mode_w_MR_h,  ac_mode_h_MR_h,  ac_mode_w_MR_c,  ac_mode_h_MR_c,  3), 
                       'BR':  make_8760_data(period_5, holiday, ac_mode_w_BR_h,  ac_mode_h_BR_h,  ac_mode_w_BR_c,  ac_mode_h_BR_c,  3), 
                       'CR1': make_8760_data(period_5, holiday, ac_mode_w_CR1_h, ac_mode_h_CR1_h, ac_mode_w_CR1_c, ac_mode_h_CR1_c, 3), 
                       'CR2': make_8760_data(period_5, holiday, ac_mode_w_CR2_h, ac_mode_h_CR2_h, ac_mode_w_CR2_c, ac_mode_h_CR2_c, 3)},
           'region6': {'MR':  make_8760_data(period_6, holiday, ac_mode_w_MR_h,  ac_mode_h_MR_h,  ac_mode_w_MR_c,  ac_mode_h_MR_c,  3), 
                       'BR':  make_8760_data(period_6, holiday, ac_mode_w_BR_h,  ac_mode_h_BR_h,  ac_mode_w_BR_c,  ac_mode_h_BR_c,  3), 
                       'CR1': make_8760_data(period_6, holiday, ac_mode_w_CR1_h, ac_mode_h_CR1_h, ac_mode_w_CR1_c, ac_mode_h_CR1_c, 3), 
                       'CR2': make_8760_data(period_6, holiday, ac_mode_w_CR2_h, ac_mode_h_CR2_h, ac_mode_w_CR2_c, ac_mode_h_CR2_c, 3)},
           'region7': {'MR':  make_8760_data(period_7, holiday, ac_mode_w_MR_h,  ac_mode_h_MR_h,  ac_mode_w_MR_c,  ac_mode_h_MR_c,  3), 
                       'BR ': make_8760_data(period_7, holiday, ac_mode_w_BR_h,  ac_mode_h_BR_h,  ac_mode_w_BR_c,  ac_mode_h_BR_c,  3), 
                       'CR1': make_8760_data(period_7, holiday, ac_mode_w_CR1_h, ac_mode_h_CR1_h, ac_mode_w_CR1_c, ac_mode_h_CR1_c, 3), 
                       'CR2': make_8760_data(period_7, holiday, ac_mode_w_CR2_h, ac_mode_h_CR2_h, ac_mode_w_CR2_c, ac_mode_h_CR2_c, 3)},
           'region8': {'MR':  make_8760_data(period_8, holiday, ac_mode_w_MR_h,  ac_mode_h_MR_h,  ac_mode_w_MR_c,  ac_mode_h_MR_c,  3), 
                       'BR ': make_8760_data(period_8, holiday, ac_mode_w_BR_h,  ac_mode_h_BR_h,  ac_mode_w_BR_c,  ac_mode_h_BR_c,  3), 
                       'CR1': make_8760_data(period_8, holiday, ac_mode_w_CR1_h, ac_mode_h_CR1_h, ac_mode_w_CR1_c, ac_mode_h_CR1_c, 3), 
                       'CR2': make_8760_data(period_8, holiday, ac_mode_w_CR2_h, ac_mode_h_CR2_h, ac_mode_w_CR2_c, ac_mode_h_CR2_c, 3)},
          }

pre_tmp = {'region1': {'MR':  make_8760_data(period_1, holiday, pre_tmp_w_MR_h,  pre_tmp_h_MR_h,  pre_tmp_w_MR_c,  pre_tmp_h_MR_c,  20.0), 
                       'BR':  make_8760_data(period_1, holiday, pre_tmp_w_BR_h,  pre_tmp_h_BR_h,  pre_tmp_w_BR_c,  pre_tmp_h_BR_c,  20.0), 
                       'CR1': make_8760_data(period_1, holiday, pre_tmp_w_CR1_h, pre_tmp_h_CR1_h, pre_tmp_w_CR1_c, pre_tmp_h_CR1_c, 20.0), 
                       'CR2': make_8760_data(period_1, holiday, pre_tmp_w_CR2_h, pre_tmp_h_CR2_h, pre_tmp_w_CR2_c, pre_tmp_h_CR2_c, 20.0)},
           'region2': {'MR':  make_8760_data(period_2, holiday, pre_tmp_w_MR_h,  pre_tmp_h_MR_h,  pre_tmp_w_MR_c,  pre_tmp_h_MR_c,  20.0), 
                       'BR':  make_8760_data(period_2, holiday, pre_tmp_w_BR_h,  pre_tmp_h_BR_h,  pre_tmp_w_BR_c,  pre_tmp_h_BR_c,  20.0), 
                       'CR1': make_8760_data(period_2, holiday, pre_tmp_w_CR1_h, pre_tmp_h_CR1_h, pre_tmp_w_CR1_c, pre_tmp_h_CR1_c, 20.0), 
                       'CR2': make_8760_data(period_2, holiday, pre_tmp_w_CR2_h, pre_tmp_h_CR2_h, pre_tmp_w_CR2_c, pre_tmp_h_CR2_c, 20.0)},
           'region3': {'MR':  make_8760_data(period_3, holiday, pre_tmp_w_MR_h,  pre_tmp_h_MR_h,  pre_tmp_w_MR_c,  pre_tmp_h_MR_c,  20.0), 
                       'BR':  make_8760_data(period_3, holiday, pre_tmp_w_BR_h,  pre_tmp_h_BR_h,  pre_tmp_w_BR_c,  pre_tmp_h_BR_c,  20.0), 
                       'CR1': make_8760_data(period_3, holiday, pre_tmp_w_CR1_h, pre_tmp_h_CR1_h, pre_tmp_w_CR1_c, pre_tmp_h_CR1_c, 20.0), 
                       'CR2': make_8760_data(period_3, holiday, pre_tmp_w_CR2_h, pre_tmp_h_CR2_h, pre_tmp_w_CR2_c, pre_tmp_h_CR2_c, 20.0)},
           'region4': {'MR':  make_8760_data(period_4, holiday, pre_tmp_w_MR_h,  pre_tmp_h_MR_h,  pre_tmp_w_MR_c,  pre_tmp_h_MR_c,  20.0), 
                       'BR':  make_8760_data(period_4, holiday, pre_tmp_w_BR_h,  pre_tmp_h_BR_h,  pre_tmp_w_BR_c,  pre_tmp_h_BR_c,  20.0), 
                       'CR1': make_8760_data(period_4, holiday, pre_tmp_w_CR1_h, pre_tmp_h_CR1_h, pre_tmp_w_CR1_c, pre_tmp_h_CR1_c, 20.0), 
                       'CR2': make_8760_data(period_4, holiday, pre_tmp_w_CR2_h, pre_tmp_h_CR2_h, pre_tmp_w_CR2_c, pre_tmp_h_CR2_c, 20.0)},
           'region5': {'MR':  make_8760_data(period_5, holiday, pre_tmp_w_MR_h,  pre_tmp_h_MR_h,  pre_tmp_w_MR_c,  pre_tmp_h_MR_c,  20.0), 
                       'BR':  make_8760_data(period_5, holiday, pre_tmp_w_BR_h,  pre_tmp_h_BR_h,  pre_tmp_w_BR_c,  pre_tmp_h_BR_c,  20.0), 
                       'CR1': make_8760_data(period_5, holiday, pre_tmp_w_CR1_h, pre_tmp_h_CR1_h, pre_tmp_w_CR1_c, pre_tmp_h_CR1_c, 20.0), 
                       'CR2': make_8760_data(period_5, holiday, pre_tmp_w_CR2_h, pre_tmp_h_CR2_h, pre_tmp_w_CR2_c, pre_tmp_h_CR2_c, 20.0)},
           'region6': {'MR':  make_8760_data(period_6, holiday, pre_tmp_w_MR_h,  pre_tmp_h_MR_h,  pre_tmp_w_MR_c,  pre_tmp_h_MR_c,  20.0), 
                       'BR':  make_8760_data(period_6, holiday, pre_tmp_w_BR_h,  pre_tmp_h_BR_h,  pre_tmp_w_BR_c,  pre_tmp_h_BR_c,  20.0), 
                       'CR1': make_8760_data(period_6, holiday, pre_tmp_w_CR1_h, pre_tmp_h_CR1_h, pre_tmp_w_CR1_c, pre_tmp_h_CR1_c, 20.0), 
                       'CR2': make_8760_data(period_6, holiday, pre_tmp_w_CR2_h, pre_tmp_h_CR2_h, pre_tmp_w_CR2_c, pre_tmp_h_CR2_c, 20.0)},
           'region7': {'MR':  make_8760_data(period_7, holiday, pre_tmp_w_MR_h,  pre_tmp_h_MR_h,  pre_tmp_w_MR_c,  pre_tmp_h_MR_c,  20.0), 
                       'BR':  make_8760_data(period_7, holiday, pre_tmp_w_BR_h,  pre_tmp_h_BR_h,  pre_tmp_w_BR_c,  pre_tmp_h_BR_c,  20.0), 
                       'CR1': make_8760_data(period_7, holiday, pre_tmp_w_CR1_h, pre_tmp_h_CR1_h, pre_tmp_w_CR1_c, pre_tmp_h_CR1_c, 20.0), 
                       'CR2': make_8760_data(period_7, holiday, pre_tmp_w_CR2_h, pre_tmp_h_CR2_h, pre_tmp_w_CR2_c, pre_tmp_h_CR2_c, 20.0)},
           'region8': {'MR':  make_8760_data(period_8, holiday, pre_tmp_w_MR_h,  pre_tmp_h_MR_h,  pre_tmp_w_MR_c,  pre_tmp_h_MR_c,  20.0), 
                       'BR':  make_8760_data(period_8, holiday, pre_tmp_w_BR_h,  pre_tmp_h_BR_h,  pre_tmp_w_BR_c,  pre_tmp_h_BR_c,  20.0), 
                       'CR1': make_8760_data(period_8, holiday, pre_tmp_w_CR1_h, pre_tmp_h_CR1_h, pre_tmp_w_CR1_c, pre_tmp_h_CR1_c, 20.0), 
                       'CR2': make_8760_data(period_8, holiday, pre_tmp_w_CR2_h, pre_tmp_h_CR2_h, pre_tmp_w_CR2_c, pre_tmp_h_CR2_c, 20.0)},
          }

vol     = {'MR': make_8760_data([1] * 365, holiday, vnt_w_MR, vnt_h_MR, [0] * 24, [0] * 24, 0.0),
           'BT': make_8760_data([1] * 365, holiday, vnt_w_BT, vnt_h_BT, [0] * 24, [0] * 24, 0.0),
           'RR': make_8760_data([1] * 365, holiday, vnt_w_RR, vnt_h_RR, [0] * 24, [0] * 24, 0.0)
          }