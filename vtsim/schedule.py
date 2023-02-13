import numpy as np

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


vnt_w_MR        = np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.250, 0.000, 0.000,  0.000, 0.000, 0.000,    
                            0.250, 0.000, 0.000,  0.000, 0.000, 0.000,  0.500, 0.500, 0.000,  0.000, 0.000, 0.000]) * 300 / 3600
vnt_h_MR        = np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.250,  0.000, 0.000, 0.000,
                            0.250, 0.000, 0.000,  0.000, 0.000, 0.500,  0.500, 0.000, 0.000,  0.000, 0.000, 0.000]) * 300 / 3600

vnt_w_BT        = np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.000,
                            0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.500, 0.250, 0.100]) * 100 / 3600
vnt_h_BT        = np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.000,
                            0.000, 0.000, 0.000,  0.000, 0.000, 0.750,  0.250, 0.000, 0.000,  0.250, 0.250, 0.100]) * 100 / 3600

vnt_w_RR        = np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.150, 0.050, 0.000,  0.020, 0.000, 0.000,
                            0.020, 0.000, 0.000,  0.000, 0.020, 0.020,  0.020, 0.020, 0.020,  0.050, 0.000, 0.070]) * 40 / 3600
vnt_h_RR        = np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.100, 0.100,  0.000, 0.030, 0.030,
                            0.000, 0.000, 0.000,  0.000, 0.050, 0.020,  0.000, 0.050, 0.020,  0.000, 0.050, 0.020]) * 40 / 3600

heat_w_MR = np.zeros(24)        #LD
heat_h_MR = np.zeros(24)
heat_w_KT = np.zeros(24)        #台所
heat_h_KT = np.zeros(24)
heat_w_BT = np.zeros(24)        #浴室
heat_h_BT = np.zeros(24)
heat_w_R1 = np.zeros(24)        #トイレ
heat_h_R1 = np.zeros(24)
heat_w_SC = np.zeros(24)        #洗面所
heat_h_SC = np.zeros(24)
heat_w_H1 = np.zeros(24)        #ホール
heat_h_H1 = np.zeros(24)
heat_w_BR = np.zeros(24)        #主寝室
heat_h_BR = np.zeros(24)
heat_w_C1 = np.zeros(24)        #子供室1
heat_h_C1 = np.zeros(24)
heat_w_C2 = np.zeros(24)        #子供室2
heat_h_C2 = np.zeros(24)

#主たる居室の発熱量#####################################################################################################
#LD人体
heat_w_MR += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  1.000, 2.000, 1.000,  1.000, 0.000, 0.000,
                       1.000, 1.000, 0.000,  0.000, 1.000, 2.000,  2.000, 3.000, 3.000,  2.000, 1.000, 1.000]) * 63.0
heat_h_MR += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 3.000,  2.000, 2.000, 2.000,
                       2.000, 1.000, 0.000,  0.000, 2.000, 3.000,  3.000, 4.000, 2.000,  2.000, 1.000, 0.000]) * 63.0
#LD照明 137.5W
heat_w_MR += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.164, 0.709, 0.382,  0.836, 0.127, 0.000,  
                       0.491, 0.382, 0.000,  0.000, 0.255, 0.509,  0.509, 0.582, 0.873,  0.509, 0.509, 0.255]) * 137.5
heat_h_MR += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.873,  1.000, 0.509, 0.509,  
                       0.745, 0.291, 0.000,  0.000, 0.509, 0.509,  0.582, 0.909, 0.509,  0.509, 0.509, 0.000]) * 137.5
#LD機器 385.08W
heat_w_MR += np.array([0.018, 0.018, 0.018,  0.018, 0.018, 0.018,  0.018, 0.543, 0.547,  0.280, 0.149, 0.018,  
                       0.280, 0.412, 0.018,  0.018, 0.280, 0.412,  0.543, 0.543, 0.543,  0.543, 0.475, 0.475]) * 385.08
heat_h_MR += np.array([0.018, 0.018, 0.018,  0.018, 0.018, 0.018,  0.018, 0.018, 0.543,  0.543, 1.000, 0.932,  
                       0.543, 0.149, 0.018,  0.018, 0.280, 0.543,  0.543, 0.280, 0.543,  0.543, 0.475, 0.018]) * 385.08

#台所の発熱量###########################################################################################################
#台所照明 36.75W
heat_w_KT += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.463, 0.463, 0.000,  0.667, 0.000, 0.000,  
                       0.925, 0.000, 0.000,  0.000, 0.925, 0.000,  0.925, 0.925, 0.925,  0.000, 0.000, 0.000]) * 36.75
heat_h_KT += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.925,  1.000, 0.000, 0.000,  
                       0.925, 0.463, 0.000,  0.000, 0.000, 0.925,  0.925, 0.925, 0.000,  0.000, 0.000, 0.000]) * 36.75
#台所機器 60W
heat_w_KT += np.array([1.000, 1.000, 1.000,  1.000, 1.000, 1.000,  1.000, 1.000, 1.000,  1.000, 1.000, 1.000, 
                       1.000, 1.000, 1.000,  1.000, 1.000, 1.000,  1.000, 1.000, 1.000,  1.000, 1.000, 1.000]) * 60.00
heat_h_KT += np.array([1.000, 1.000, 1.000,  1.000, 1.000, 1.000,  1.000, 1.000, 1.000,  1.000, 1.000, 1.000, 
                       1.000, 1.000, 1.000,  1.000, 1.000, 1.000,  1.000, 1.000, 1.000,  1.000, 1.000, 1.000]) * 60.00
#台所機器 34.76W
heat_w_KT += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.500, 0.000, 0.000,  0.000, 0.000, 0.000, 
                       0.500, 0.000, 0.000,  0.000, 0.000, 0.000,  1.000, 0.000, 0.000,  0.000, 0.000, 0.000]) * 34.76
heat_h_KT += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.500,  0.000, 0.000, 0.000, 
                       0.500, 0.000, 0.000,  0.000, 0.000, 1.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.000]) * 34.76

#浴室の発熱量###########################################################################################################
#浴室照明 40.5W
heat_w_BT += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  
                       0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.333,  0.667, 1.000, 0.000]) * 40.5
heat_h_BT += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  
                       0.000, 0.000, 0.000,  0.000, 0.667, 0.333,  0.000, 0.000, 0.000,  0.667, 1.000, 0.000]) * 40.5

#トイレの発熱量#########################################################################################################              
#トイレ機器 30.0W
heat_w_R1 += np.array([1.000, 1.000, 1.000,  1.000, 1.000, 1.000,  1.000, 1.000, 1.000,  1.000, 1.000, 1.000,  
                       1.000, 1.000, 1.000,  1.000, 1.000, 1.000,  1.000, 1.000, 1.000,  1.000, 1.000, 1.000]) * 30.0
heat_h_R1 += np.array([1.000, 1.000, 1.000,  1.000, 1.000, 1.000,  1.000, 1.000, 1.000,  1.000, 1.000, 1.000,  
                       1.000, 1.000, 1.000,  1.000, 1.000, 1.000,  1.000, 1.000, 1.000,  1.000, 1.000, 1.000]) * 30.0
#トイレ照明 8.55W
heat_w_R1 += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  1.000, 0.333, 0.000,  0.111, 0.000, 0.000,  
                       0.111, 0.000, 0.000,  0.000, 0.111, 0.111,  0.111, 0.111, 0.111,  0.333, 0.000, 0.444]) * 8.55
heat_h_R1 += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.667, 0.667,  0.000, 0.222, 0.222,  
                       0.000, 0.000, 0.000,  0.000, 0.333, 0.111,  0.000, 0.333, 0.111,  0.000, 0.333, 0.111]) * 8.55

#洗面室の発熱###########################################################################################################
#洗面脱衣照明 66.5W
heat_w_SC += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.286, 0.571, 0.238,  0.524, 0.286, 0.000,  
                       0.000, 0.286, 0.000,  0.000, 0.095, 0.095,  0.190, 0.286, 0.214,  1.000, 0.929, 0.286]) * 66.5
heat_h_SC += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.786, 0.786,  0.857, 0.000, 0.095,  
                       0.000, 0.000, 0.000,  0.000, 0.452, 0.500,  0.190, 0.000, 0.000,  0.714, 0.929, 0.286]) * 66.5
#洗面脱衣機器 118.75W
heat_w_SC += np.array([0.097, 0.097, 0.097,  0.097, 0.097, 0.097,  0.097, 0.548, 0.227,  0.097, 0.097, 0.097,  
                       0.097, 0.097, 0.097,  0.097, 0.097, 0.097,  0.097, 0.097, 0.097,  1.000, 0.097, 0.548]) * 118.75
heat_h_SC += np.array([0.097, 0.097, 0.097,  0.097, 0.097, 0.097,  0.097, 0.548, 0.227,  0.097, 0.097, 0.097,  
                       0.097, 0.097, 0.097,  0.097, 0.097, 0.548,  0.097, 0.097, 0.097,  0.548, 0.097, 0.548]) * 118.75

#ホールの発熱###########################################################################################################
#ホール照明 57W
heat_w_H1 += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.500, 1.000, 1.000,  1.000, 0.500, 0.000,
                       0.000, 0.000, 0.000,  0.000, 1.000, 1.000,  1.000, 1.000, 1.000,  1.000, 1.000, 0.500]) * 57.0
heat_h_H1 += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.750, 1.000,  1.000, 1.000, 1.000,  
                       1.000, 0.250, 0.000,  0.000, 0.000, 0.000,  0.500, 1.000, 1.000,  1.000, 1.000, 0.250]) * 57.0
#ホール照明 114W
heat_w_H1 += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.250, 0.500, 0.250,  0.500, 0.250, 0.000,
                       0.000, 0.250, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.250,  1.000, 1.000, 0.250]) * 114.0
heat_h_H1 += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.750, 0.750,  1.000, 0.000, 0.000,  
                       0.000, 0.000, 0.000,  0.000, 0.500, 0.250,  0.000, 0.000, 0.000,  0.250, 0.250, 0.250]) * 114.0

#寝室の発熱量###########################################################################################################
#寝室人体
heat_w_BR += np.array([2.000, 2.000, 2.000,  2.000, 2.000, 2.000,  1.000, 0.000, 0.000,  0.000, 0.000, 0.000,
                       0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 1.000]) * 63.0
heat_h_BR += np.array([2.000, 2.000, 2.000,  2.000, 2.000, 2.000,  2.000, 1.000, 0.000,  0.000, 0.000, 0.000,
                       0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 2.000]) * 63.0
#寝室照明 52.5W
heat_w_BR += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.667, 0.000, 0.000,  0.000, 0.000, 0.000,
                       0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.000]) * 52.5
heat_h_BR += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  1.000, 0.000, 0.000,  0.000, 0.000, 0.000,
                       0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.000]) * 52.5
#寝室照明 412.5W
heat_w_BR += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.667, 0.000, 0.000,  0.000, 0.000, 0.000,
                       0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.000]) * 412.5
heat_h_BR += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  1.000, 0.000, 0.000,  0.000, 0.000, 0.000,
                       0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.000]) * 412.5

#子供室1の発熱量########################################################################################################
#子供室1人体
heat_w_BR += np.array([1.000, 1.000, 1.000,  1.000, 1.000, 1.000,  1.000, 0.000, 0.000,  0.000, 0.000, 0.000,
                       0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 1.000,  0.000, 1.000, 1.000]) * 63.0
heat_h_BR += np.array([1.000, 1.000, 1.000,  1.000, 1.000, 1.000,  1.000, 1.000, 1.000,  1.000, 1.000, 1.000,
                       0.000, 0.000, 0.000,  0.000, 1.000, 1.000,  1.000, 0.000, 1.000,  1.000, 1.000, 1.000]) * 63.0
#子供室1照明 70W
heat_w_BR += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.500, 0.000, 0.000,
                       0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.750,  0.250, 1.000, 1.000]) * 70.0
heat_h_BR += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.750, 1.000, 1.000,
                       0.000, 0.000, 0.000,  0.000, 1.000, 1.000,  0.500, 0.000, 1.000,  0.250, 1.000, 0.000]) * 70.0
#子供室1照明 80W
heat_w_BR += np.array([0.188, 0.188, 0.188,  0.188, 0.188, 0.188,  0.188, 0.188, 0.188,  0.188, 0.188, 0.188,
                       0.188, 0.188, 0.188,  0.188, 0.188, 0.188,  0.188, 0.188, 0.750,  0.375, 1.000, 0.438]) * 80.0
heat_h_BR += np.array([0.188, 0.188, 0.188,  0.188, 0.188, 0.188,  0.188, 0.188, 0.188,  0.797, 1.000, 1.000,
                       0.188, 0.188, 0.188,  0.188, 0.250, 0.250,  0.219, 0.188, 1.000,  0.391, 1.000, 0.188]) * 80.0

#子供室2の発熱量########################################################################################################
#子供室2人体
heat_w_BR += np.array([1.000, 1.000, 1.000,  1.000, 1.000, 1.000,  1.000, 0.000, 0.000,  0.000, 0.000, 0.000,
                       0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  1.000, 0.000, 0.000,  1.000, 1.000, 1.000]) * 63.0
heat_h_BR += np.array([1.000, 1.000, 1.000,  1.000, 1.000, 1.000,  1.000, 1.000, 0.000,  1.000, 1.000, 1.000,
                       1.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 1.000,  1.000, 1.000, 1.000]) * 63.0
#子供室2照明 70W
heat_w_BR += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.500, 0.000, 0.000,
                       0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.500, 0.500, 0.000,  0.750, 1.000, 0.250]) * 70.0
heat_h_BR += np.array([0.000, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 0.250,  1.000, 1.000, 1.000,
                       0.500, 0.000, 0.000,  0.000, 0.000, 0.000,  0.000, 0.000, 1.000,  1.000, 1.000, 0.000]) * 70.0
#子供室2照明 50W
heat_w_BR += np.array([0.060, 0.060, 0.060,  0.060, 0.060, 0.060,  0.060, 0.060, 0.060,  0.060, 0.060, 0.060,
                       0.060, 0.060, 0.060,  0.060, 0.060, 0.060,  0.060, 0.060, 0.060,  0.000, 0.765, 0.295]) * 50.0
heat_h_BR += np.array([0.060, 0.060, 0.060,  0.060, 0.060, 0.060,  0.060, 0.060, 0.060,  0.060, 0.060, 0.060,
                       0.060, 0.060, 0.060,  0.060, 0.060, 0.060,  0.060, 0.060, 1.000,  1.000, 1.000, 0.060]) * 50.0


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

heat    ={'LD':            make_8760_data([1] * 365, holiday, heat_w_MR.tolist(), heat_w_MR.tolist(), [0] * 24, [0] * 24, 0.0),
          'Kitchen':       make_8760_data([1] * 365, holiday, heat_w_KT.tolist(), heat_w_KT.tolist(), [0] * 24, [0] * 24, 0.0),
          'Bath':          make_8760_data([1] * 365, holiday, heat_w_BT.tolist(), heat_w_BT.tolist(), [0] * 24, [0] * 24, 0.0),
          'Toilet':        make_8760_data([1] * 365, holiday, heat_w_R1.tolist(), heat_w_R1.tolist(), [0] * 24, [0] * 24, 0.0),
          'Sanitary':      make_8760_data([1] * 365, holiday, heat_w_SC.tolist(), heat_w_SC.tolist(), [0] * 24, [0] * 24, 0.0),
          'Hall':          make_8760_data([1] * 365, holiday, heat_w_H1.tolist(), heat_w_H1.tolist(), [0] * 24, [0] * 24, 0.0),
          'Bedroom':       make_8760_data([1] * 365, holiday, heat_w_BR.tolist(), heat_w_BR.tolist(), [0] * 24, [0] * 24, 0.0),
          'ChildrenRoom1': make_8760_data([1] * 365, holiday, heat_w_C1.tolist(), heat_w_C1.tolist(), [0] * 24, [0] * 24, 0.0),
          'ChildrenRoom2': make_8760_data([1] * 365, holiday, heat_w_C2.tolist(), heat_w_C2.tolist(), [0] * 24, [0] * 24, 0.0)}
