class fan_spec:
    def __init__(self):
        self.qmax = {}
        self.pmax = {}
        self.q1   = {}
        self.p1   = {}

    def set_spec(self, name, qmax, pmax, q1, p1):
        self.qmax[name] = qmax
        self.pmax[name] = pmax
        self.q1[name]   = q1
        self.p1[name]   = p1

    def fan_dict(self, df, dic):
        qmax_d, pmax_d, q1_d, p1_d = {}, {}, {}, {}

        for k in dic:  
            qmax_d[k] = self.qmax[dic[k]]
            pmax_d[k] = self.pmax[dic[k]]
            q1_d[k]   = self.q1[dic[k]]
            p1_d[k]   = self.p1[dic[k]]

        qmax = df.replace(qmax_d).to_list()
        pmax = df.replace(pmax_d).to_list()
        q1   = df.replace(q1_d).to_list()
        p1   = df.replace(p1_d).to_list()
        
        return {'qmax': qmax, 'pmax': pmax, 'q1': q1, 'p1': p1}

fs  = fan_spec()

fs.set_spec('YCC_L',  100 / 3600, 110.0, 100 / 3600,  95.0)
fs.set_spec('YCC_M',  200 / 3600, 110.0, 200 / 3600,  80.0)
fs.set_spec('YCC_H',  390 / 3600, 110.0, 200 / 3600,  80.0)

fs.set_spec('RARA_mini', 70 / 3600,  10.0,   0 / 3600, 10.0)
fs.set_spec('RARA0',    100 / 3600,  10.0,   0 / 3600, 10.0)
fs.set_spec('RARA1',    150 / 3600,  80.0,   0 / 3600, 80.0)
fs.set_spec('RARA2',    210 / 3600, 120.0,   0 / 3600, 120.0)

fs.set_spec('SUMIKA', 200 / 3600, 100.0, 200 / 3600, 100.0)

fs.set_spec('AEROTHECH_50Hz', 1850 / 3600, 440.0, 1640 / 3600, 100.0)
fs.set_spec('AEROTHECH_40Hz', 1550 / 3600, 340.0, 1290 / 3600, 100.0)
fs.set_spec('AEROTHECH_25Hz', 1150 / 3600, 300.0,  840 / 3600, 100.0)
fs.set_spec('AEROTHECH_15Hz',  680 / 3600, 180.0,    0 / 3600, 180.0)
fs.set_spec('AEROTHECH_10Hz',  400 / 3600,  80.0,    0 / 3600,  80.0)

#旧バージョン
#YCC Fan Spec
qmax_YCC_L, pmax_YCC_L, q1_YCC_L, p1_YCC_L = 100 / 3600, 110.0, 100 / 3600,  95.0
qmax_YCC_M, pmax_YCC_M, q1_YCC_M, p1_YCC_M = 200 / 3600, 110.0, 200 / 3600,  80.0
qmax_YCC_H, pmax_YCC_H, q1_YCC_H, p1_YCC_H = 390 / 3600, 110.0, 200 / 3600,  80.0

#RARA Fan Spec
qmax_RARA0, pmax_RARA0, q1_RARA0, p1_RARA0 = 100 / 3600,  10.0,   0 / 3600,  10.0
qmax_RARA1, pmax_RARA1, q1_RARA1, p1_RARA1 = 150 / 3600,  80.0,   0 / 3600,  80.0
qmax_RARA2, pmax_RARA2, q1_RARA2, p1_RARA2 = 210 / 3600, 120.0,   0 / 3600, 120.0