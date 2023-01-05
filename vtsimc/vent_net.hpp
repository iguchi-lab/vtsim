using namespace std;

class Vent_Net{
public:
    int i, i1, i2, vn_type;                                         //入口出口のノード番号、換気回路網の種類
    vector<double> h1, h2;                                          //入口出口の高さ
    vector<double> alpha, area;                                     //単純開口　開口率、面積
    vector<double> a, n;                                            //隙間　　　隙間量、隙間特性値
    vector<double> qv, qt;                                          //風量、移流に伴う熱量
    vector<double> eta;                                             //粉塵除去率
    vector<double> q_max, p_max, q1, p1;                            //送風ファンの風量・圧力

    Vent_Net(long length, int i_, int i1_, int i2_, int vn_type_, vector<double> h1_, vector<double> h2_){
        i       = i_;
        i1      = i1_;
        i2      = i2_;
        vn_type = vn_type_;
        h1      = h1_;
        h2      = h2_; 
        qv.assign(length, 0.0);                                     //風量を0に初期化
        qt.assign(length, 0.0);                                     //熱量を0に初期化
    }

    double get_qv(double dp, long ts){
        switch(vn_type){
            case VN_SIMPLE:                                                         //単純開口
                if(dp >= 0)                                   return  alpha[ts] * area[ts] * sqrt( 2.0 * dp / RHO20);
                else                                          return -alpha[ts] * area[ts] * sqrt(-2.0 * dp / RHO20);
            case VN_GAP:                                                            //隙間
                if(dp >= 0)                                   return  a[ts] * pow( dp, 1 / n[ts]);
                else                                          return -a[ts] * pow(-dp, 1 / n[ts]);
            case VN_FIX:
            case VN_AIRCON:                                   return qv[ts];      //固定・エアコン 
            case VN_FAN:                                                          //送風ファン
                if (-dp <= 0.0)                               return q_max[ts];
                else if((0 < -dp) && (-dp <= p1[ts]))         return q_max[ts] + (q1[ts] - q_max[ts]) / p1[ts] * (-dp);
                else if((p1[ts] < -dp) && (-dp <= p_max[ts])) return q1[ts] - q1[ts] * (p_max[ts] - p1[ts]) * (-dp - q1[ts]);
                else                                          return 0.0;
            default:                                          return 0.0;
        }
    }

    double get_qt(double dt, long ts){
        switch(vn_type){
            case VN_AIRCON:                                 return 0.0;           //エアコンの場合のみ0
            default:                                        return RHO20 * AIR_CP * qv[ts] * dt;
        }
    }
};