import lib.mat_data
import lib.wal_data

class Material:
    def set_material(self, name, lambda_, capa_):
        self.material[name] = {'lambda': lambda_, 'capa': capa_}

    def __init__(self):
        self.material = mat_data.data

class Wall:
    def __init__(self, wall):
        
        r = [0.0] * len(wall)
        u = [0.0] * len(wall)
        c = [0.0] * len(wall)
        
        mat = Material()

        for i, l in enumerate(wall):
            for w in l[0]:
                if w[0]  == '通気層': 
                    r[i] += 0.110
                    c[i] += w[1] * 1.006 * 1.2 * 1000            #m  *  kJ / (kg・K)  *  kg / m3  *  J / kJ
                elif w[0] == '中空層':
                    r[i] += 0.090
                    c[i] += w[1] * 1.006 * 1.2 * 1000            #m  *  kJ / (kg・K)  *  ka / m3  *  J / kJ
                else:
                    r[i] += w[1] / mat.material[w[0]]['lambda']
                    c[i] += w[1] * mat.material[w[0]]['capa'] * 1000 #J/(m2・K)

            u[i] = 1 / r[i] * l[1]
            c[i] = c[i] * l[1]

        self.u_value = sum(u)
        self.capa_w  = sum(c)
    
    def spec(self):
        return {'U_w': self.u_value, 'capa_w': self.capa_w, 'eta_w': 0.8}

wall = {}
for wl in wal_data.basic:   wall[wl] = Wall(wal_data.basic[wl])
for wl in wal_data.kameido: wall[wl] = Wall(wal_data.kameido[wl])
for wl in wal_data.FPJ:     wall[wl] = Wall(wal_data.FPJ[wl])
for wl in wal_data.okayama: wall[wl] = Wall(wal_data.okayama[wl])
for wl in wal_data.oosaka:  wall[wl] = Wall(wal_data.oosaka[wl])