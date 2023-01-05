#include "vtsim_set.hpp"            //計算の設定　ヘッダーファイルの読み込み
#include "node.hpp"                 //ノード　ヘッダーファイルの読み込み
#include "vent_net.hpp"             //換気回路網　ヘッダーファイルの読み込み
#include "thrm_net.hpp"             //熱回路網　ヘッダーファイルの読み込み
#include <iostream>
#include <fstream>

//#define DEBUG_ON

#ifdef  DEBUG_ON
#define LOG_PRINT(...)     ofs << __FILE__ << " (" << __LINE__ << ") " << __func__ << ":" << __VA_ARGS__
#define LOG_CONTENTS(...)  ofs << __VA_ARGS__
const char *logfileName = "log.txt";
ofstream ofs(logfileName);
#else
#define LOG_PRINT(...)
#define LOG_CONTENTS(...)
#endif

#include "mymath.hpp"

using namespace::std;

double calc_x_hum(double t, double h){
    double f  = - 7.90295 * (373.16 / (t + 273.15) - 1) + 
                  5.02808 * log10(373.16 / (t + 273.15)) - 
                  1.3816e-7 * (pow(10, 11.344  * (1 - (t + 273.15) / 373.16)) - 1) + 
                  8.1328e-3 * (pow(10, -3.4919 * (373.16 / (t + 273.15) - 1)) - 1) + log10(1013.246);
    double ps = pow(10, f) * 100;                                              //飽和水蒸気圧
    double e  = h / 100 * ps;                                                    //水蒸気圧
    double x  = 0.622 * e / (101325 - e);                                        //絶対湿度
    
    return x;
}

class CalcStatus{
public:
    long length      = 0;
    double t_step    = 3600;
    int solve        = SOLVE_LU;
    double step_p    = STEP_P, 
           vent_err  = VENT_ERR, 
           step_t    = STEP_T, 
           thrm_err  = THRM_ERR, 
           conv_err  = CONV_ERR, 
           sor_ratio = SOR_RATIO, 
           sor_err   = SOR_ERR;
};

class VTSim{
public:
    CalcStatus sts;
    vector<Node> sn;                                        //ノード
    map<string, int>  node, v_net, t_net;
    vector<Vent_Net>  vn;                                   //換気回路網
    vector<Thrm_Net>  tn;                                   //熱回路網
    vector<int> v_idc, c_idc, t_idc, x_idc;
    vector<int> i_vn_ac, i_tn_ac;

    void setup(CalcStatus sts_){
        sts = sts_;
        sn.clear();
        node.clear();
        v_net.clear();
        t_net.clear();
        vn.clear();
        tn.clear();
        v_idc.clear();
        c_idc.clear();
        t_idc.clear();
        x_idc.clear();
        i_vn_ac.clear();
        i_tn_ac.clear();
    }

    void sn_add(int i, tuple<int, int, int, int> flag){
        sn.push_back(Node(sts.length, i, flag));
        if(get<0>(flag) == SN_CALC)     v_idc.push_back(i);
        if(get<1>(flag) == SN_CALC)     c_idc.push_back(i);
        if(get<2>(flag) == SN_CALC)     t_idc.push_back(i);
        if(get<3>(flag) == SN_CALC)     x_idc.push_back(i);
    }

    void set_node(string name, int i){
        node[name] = i;
    }

    void vn_add(int i, string name, int i1, int i2, int vn_type, vector<double> h1, vector<double> h2){
        v_net[name] = i;
        vn.push_back(Vent_Net(sts.length, i, i1, i2, vn_type, h1, h2));
    }

    void tn_add(int i, string name, int i1, int i2, int tn_type){
        t_net[name] = i;
        tn.push_back(Thrm_Net(sts.length, i, i1, i2, tn_type));
    }
    
    void vn_aircon_add(int i){
        i_vn_ac.push_back(i);
    }

    void tn_aircon_add(int i){
        i_tn_ac.push_back(i);
    }
    
    void change_sn_t_flag(int i, int flag_){
        get<2>(sn[i].flag) = flag_;
        t_idc.clear();
        for(unsigned int i = 0; i < sn.size(); i++)   
            if(get<2>(sn[i].flag) == SN_CALC)     t_idc.push_back(sn[i].i);
    }

    void set_aircon_status(long ts, int i, int flag_){
        if(flag_ == AC_ON){
            tn[i_tn_ac[i]].aircon_on[ts] = AC_ON;
            change_sn_t_flag(tn[i_tn_ac[i]].i1, SN_FIX);
            vn[i_vn_ac[i]].vn_type = VN_AIRCON;
        }
        else if(flag_ == AC_OFF){
            tn[i_tn_ac[i]].aircon_on[ts] = AC_OFF;
            change_sn_t_flag(tn[i_tn_ac[i]].i1, SN_CALC);
            vn[i_vn_ac[i]].vn_type = VN_FIX;
        }
    }

    vector<double> qv_sum(vector<double> p, long ts){
        vector<double> qvsum(sn.size(), 0.0);                                                                             //風量収支の初期化                                 
        for(unsigned int i = 0; i < vn.size(); i++){
            double rgh1 = get_rho(sn[vn[i].i1].t[ts]) * G * vn[i].h1[ts];
            double rgh2 = get_rho(sn[vn[i].i2].t[ts]) * G * vn[i].h2[ts];
            double qv   = vn[i].get_qv((p[vn[i].i1] - rgh1) - (p[vn[i].i2] - rgh2), ts);
            qvsum[vn[i].i1] -= qv;                                                                                        //風量収支の計算
            qvsum[vn[i].i2] += qv;                                                                                        //風量収支の計算      
        }
        return qvsum;
    }

    double calc_vent(long ts){
        vector<double>          p0(sn.size()), qvsum_0(sn.size()), qvsum_d(sn.size());
        vector<vector<double>>  a(v_idc.size(), vector<double>(v_idc.size()));
        vector<double>          b(v_idc.size()), dp(v_idc.size());
        double                  rmse;
        int                     itr = 0;
        
        for(unsigned int i = 0; i < sn.size(); i++)    p0[sn[i].i] = sn[i].p[ts];                                                //圧力の初期化       
        do{
            qvsum_0 = qv_sum(p0, ts);

            for(unsigned int j = 0; j < v_idc.size(); j++){
                p0[v_idc[j]] += sts.step_p;                                                                             //ダミー圧力の作成       
                qvsum_d = qv_sum(p0, ts);                                                                            //ダ三－風量収支の計算
                for(unsigned int i = 0; i < v_idc.size(); i++)  a[i][j] = (qvsum_d[v_idc[i]] - qvsum_0[v_idc[i]]) / sts.step_p;  //aの計算
                b[j] = -qvsum_0[v_idc[j]];
                p0[v_idc[j]] -= sts.step_p;
            }

            if(sts.solve == SOLVE_SOR)  dp = SOR(a, b, v_idc.size(), sts.sor_ratio, sts.sor_err);                       //SOR法による計算     
            else                        dp = LU(a, b, v_idc.size()); 

            for(unsigned int i = 0; i < v_idc.size(); i++)   p0[v_idc[i]] += dp[i];                                     //圧力の更新

            qvsum_0 = qv_sum(p0, ts);
            rmse = 0.0;    
            for(unsigned int i = 0; i < v_idc.size(); i++)  rmse += pow(qvsum_0[v_idc[i]], 2.0) / v_idc.size();
            LOG_PRINT(itr << ": rmse = " << sqrt(rmse) << endl);
            itr++;

        }while(sts.vent_err < sqrt(rmse));

        for(unsigned int i = 0; i < v_idc.size(); i++)   sn[v_idc[i]].p[ts] = p0[v_idc[i]];                             //圧力の計算                                                               
        return sqrt(rmse);
    }

    vector<double> qt_sum(vector<double> t, long ts){
        vector<double> qtsum(sn.size(), 0.0);           

        for(unsigned int i = 0; i < vn.size(); i++){                                                                    //移流に伴う熱移動
            if(vn[i].qv[ts] > 0)            qtsum[vn[i].i2] += vn[i].get_qt(t[vn[i].i1] - t[vn[i].i2], ts);
            else                            qtsum[vn[i].i1] += vn[i].get_qt(t[vn[i].i1] - t[vn[i].i2], ts);
        }

        for(unsigned int i = 0; i < tn.size(); i++){                                                                    //貫流、日射、発熱による熱移動
            switch(tn[i].tn_type){
                case TN_SIMPLE:
                case TN_GROUND:
                    qtsum[tn[i].i1] -= tn[i].get_qt(t[tn[i].i1] - t[tn[i].i2], ts);     
                    qtsum[tn[i].i2] += tn[i].get_qt(t[tn[i].i1] - t[tn[i].i2], ts);
                    break;
                case TN_SOLAR:
                    qtsum[tn[i].i1] += tn[i].get_qt_s(sn[tn[i].i2].h_sr[ts], ts);
                    break;
                case TN_HEATER:
                    qtsum[tn[i].i1] += tn[i].get_qt_h(sn[tn[i].i2].h_inp[ts], ts);
                    break;
            }    
        }

        for(unsigned int i = 0; i < i_tn_ac.size(); i++){
            if(tn[i_tn_ac[i]].aircon_on[ts] == AC_ON){                                                                  //エアコンによる温度調整
                qtsum[tn[i_tn_ac[i]].i2] -= qtsum[tn[i_tn_ac[i]].i1];
                qtsum[tn[i_tn_ac[i]].i1]  = 0.0;
            }
        }

        return qtsum;
    }

    double calc_thrm(long ts){
        vector<double>          t0(sn.size()), qtsum_0(sn.size()), qtsum_d(sn.size());
        vector<vector<double>>  a(t_idc.size() + i_tn_ac.size(), vector<double>(t_idc.size() + i_tn_ac.size()));                                  //エアコン分も余計に確保
        vector<double>          b(t_idc.size() + i_tn_ac.size()), dt(t_idc.size() + i_tn_ac.size());                                              //エアコン分も余計に確保
        double                  rmse;        
        int                     itr = 0;

        for(unsigned int i = 0; i < sn.size(); i++)     t0[sn[i].i] = sn[i].t[ts];                                      //温度の初期化

        do{
            for(unsigned int i = 0; i < i_tn_ac.size(); i++){
                if(itr == 0)    
                    set_aircon_status(ts, i, AC_OFF);   //初回は暖冷房OFFで計算し、暖冷房が必要かどうかチェックする
                else if(((tn[i_tn_ac[i]].ac_mode[ts] == AC_AUTO)) ||
                        ((tn[i_tn_ac[i]].ac_mode[ts] == AC_HEATING) && (t0[tn[i_tn_ac[i]].i1] <= tn[i_tn_ac[i]].pre_tmp[ts])) ||
                        ((tn[i_tn_ac[i]].ac_mode[ts] == AC_COOLING) && (t0[tn[i_tn_ac[i]].i1] >= tn[i_tn_ac[i]].pre_tmp[ts]))){
                    set_aircon_status(ts, i, AC_ON);    //設定温度を満たしていない場合は、暖冷房ON
                }
            }
            
            qtsum_0 = qt_sum(t0, ts);                                                                                   //熱量収支の計算

            for(unsigned int j = 0; j < t_idc.size(); j++){
                t0[t_idc[j]] += sts.step_t;                                                                             //ダミー温度の作成
                qtsum_d = qt_sum(t0, ts);                                                                               //ダミー熱量の計算
                for(unsigned int i = 0; i < t_idc.size(); i++)  a[i][j] = (qtsum_d[t_idc[i]] - qtsum_0[t_idc[i]]) / sts.step_t;  //aの計算
                b[j] = -qtsum_0[t_idc[j]];                                                                              //bの計算
                t0[t_idc[j]] -= sts.step_t;                                                                             //ダミー温度を戻す
            }

            if(sts.solve == SOLVE_SOR)  dt = SOR(a, b, t_idc.size(), sts.sor_ratio, sts.sor_err);                       //SOR法による計算
            else                        dt = LU(a, b, t_idc.size());

            for(unsigned int i = 0; i < t_idc.size(); i++)   t0[t_idc[i]] += dt[i];                                     //温度の更新

            for(unsigned int i = 0; i < i_tn_ac.size(); i++){
                double  pre_t = tn[i_tn_ac[i]].pre_tmp[ts], t_in = t0[tn[i_tn_ac[i]].i1], t_out = t0[tn[i_tn_ac[i]].i2];

                if(((tn[i_tn_ac[i]].ac_mode[ts] == AC_AUTO)    && (abs(t_in - pre_t) > sts.thrm_err)) ||
                   ((tn[i_tn_ac[i]].ac_mode[ts] == AC_HEATING) && (    t_in          < pre_t       )) ||
                   ((tn[i_tn_ac[i]].ac_mode[ts] == AC_COOLING) && (    t_in          > pre_t       ))){     //設定温度を満たしているかチェック
                    t0[tn[i_tn_ac[i]].i1] = pre_t;
                }
                else{       //設定温度を満たしていれば最大能力を超えているかチェック & 除湿による吹き出し口絶対湿度の計算
                    if(tn[i_tn_ac[i]].ac_model == AC_RAC){
                        double  x_in  = sn[tn[i_tn_ac[i]].i1].x[ts], x_out = sn[tn[i_tn_ac[i]].i2].x[ts];
                        double  qv    = vn[i_vn_ac[i]].qv[ts], To = tn[i_tn_ac[i]].To[ts], ho = tn[i_tn_ac[i]].ho[ts];
                        double  t_out_lim, x_out_lim;
                        if(t_in < t_out){
                            sn[tn[i_tn_ac[i]].i2].x[ts] = x_in;
                            LOG_PRINT("RAC:" << pre_t << ", " << t_in << ", " <<  t_out << ", " << 
                                                x_in << ", " << x_out << ", " << qv << ", " << To << ", " << ho << endl);
                            tie(t_out_lim, x_out_lim) = tn[i_tn_ac[i]].get_RAC_hs_out(t_in, t_out, x_in, x_out, qv, To, ho, ts);
                            if(t_out - t_out_lim > 0.001){
                                tn[i_tn_ac[i]].pre_tmp[ts] -= 0.01;
                                t0[tn[i_tn_ac[i]].i1]      -= 0.01;
                            }
                        }
                        if(t_in > t_out){
                            sn[tn[i_tn_ac[i]].i2].x[ts] = min(x_out, calc_x_hum(t_out, 95.0));
                            LOG_PRINT("RAC:" << pre_t << ", " << t_in << ", " <<  t_out << ", " << 
                                       x_in << ", " << min(x_out, calc_x_hum(t_out, 95.0)) << ", " << qv << ", " << To << ", " << ho << endl);
                            tie(t_out_lim, x_out_lim) = tn[i_tn_ac[i]].get_RAC_hs_out(t_in, t_out, x_in, min(x_out, calc_x_hum(t_out, 95.0)), qv, To, ho, ts);
                            if(t_out_lim - t_out > 0.001 || x_out_lim - x_out > 0.0001){
                                tn[i_tn_ac[i]].pre_tmp[ts]  += 0.01;
                                t0[tn[i_tn_ac[i]].i1]       += 0.01;
                            }
                        }
                    }
                    if(tn[i_tn_ac[i]].ac_model == AC_DUCT_C){ 
                        double  x_in  = sn[tn[i_tn_ac[i]].i1].x[ts], x_out = sn[tn[i_tn_ac[i]].i2].x[ts];
                        double  qv    = vn[i_vn_ac[i]].qv[ts], To = tn[i_tn_ac[i]].To[ts], ho = tn[i_tn_ac[i]].ho[ts];
                        double  t_out_lim, x_out_lim;      
                        if(t_in < t_out){
                            sn[tn[i_tn_ac[i]].i2].x[ts] = x_in;
                            LOG_PRINT("RAC:" << pre_t << ", " << t_in << ", " <<  t_out << ", " << 
                                                x_in << ", " << x_out << ", " << qv << ", " << To << ", " << ho << endl);
                            tie(t_out_lim, x_out_lim) 
                                = tn[i_tn_ac[i]].get_DUCT_C_hs_out(t_in, t_out, x_in, x_out, qv, To, ho, ts);
                        }
                        if(t_in > t_out){
                            sn[tn[i_tn_ac[i]].i2].x[ts] = min(x_out, calc_x_hum(t_out, 95.0));
                            LOG_PRINT("DUCT_C:" << pre_t << ", " << t_in << ", " <<  t_out << ", " << 
                                       x_in << ", " << min(x_out, calc_x_hum(t_out, 95.0)) << ", " << qv << ", " << To << ", " << ho << endl);
                            tie(t_out_lim, x_out_lim) = tn[i_tn_ac[i]].get_DUCT_C_hs_out(t_in, t_out, x_in, min(x_out, calc_x_hum(t_out, 95.0)), qv, To, ho, ts);
                        }
                    }
                }
            }

            qtsum_0 = qt_sum(t0, ts);
            rmse = 0.0;
            for(unsigned int i = 0; i < t_idc.size(); i++)   rmse += pow(qtsum_0[t_idc[i]], 2.0) / t_idc.size(); 
            LOG_PRINT(itr << ": rmse = " << sqrt(rmse) << endl);
            itr++;

        }while(sts.thrm_err < sqrt(rmse));

        for(unsigned int i = 0; i < t_idc.size(); i++)      sn[t_idc[i]].t[ts] = t0[t_idc[i]];                          //温度の計算
        for(unsigned int i = 0; i < i_tn_ac.size(); i++){    
            sn[tn[i_tn_ac[i]].i1].t[ts] = t0[tn[i_tn_ac[i]].i1];        //エアコン設定温度ノードの温度を更新
            sn[tn[i_tn_ac[i]].i2].t[ts] = t0[tn[i_tn_ac[i]].i2];        //エアコン設定温度ノードの温度を更新
        }
        return sqrt(rmse);
    }

    void calc_qv(long ts1, long ts2){
        //LOG_PRINT("ts1 = " << ts1 << " <<<<---- " << "ts2 = " <<  ts2 << endl);
        for(unsigned int i = 0; i < vn.size(); i++){    
            double rgh1 = get_rho(sn[vn[i].i1].t[ts2]) * G * vn[i].h1[ts2];
            double rgh2 = get_rho(sn[vn[i].i2].t[ts2]) * G * vn[i].h2[ts2];
            if(vn[i].vn_type != VN_FIX)
                vn[i].qv[ts1] = vn[i].get_qv((sn[vn[i].i1].p[ts2] - rgh1) - (sn[vn[i].i2].p[ts2] - rgh2), ts2);             //風量の計算
        }
    }

    void calc_qt(long ts1, long ts2){
        //LOG_PRINT("ts1 = " << ts1 << " <<<<---- " << "ts2 = " << ts2 << endl);
        for(unsigned int i = 0; i < tn.size(); i++){    
            switch(tn[i].tn_type){
                case TN_SIMPLE: 
                case TN_GROUND:         tn[i].qt[ts1] =   tn[i].get_qt(sn[tn[i].i1].t[ts2] - sn[tn[i].i2].t[ts2], ts2);
                                        break;
                case TN_SOLAR:          tn[i].qt[ts1] = - tn[i].get_qt_s(sn[tn[i].i2].h_sr[ts2], ts2);
                                        break;
                case TN_HEATER:         tn[i].qt[ts1] = - tn[i].get_qt_h(sn[tn[i].i2].h_inp[ts2], ts2);
                                        break;
            }
        }
        for(unsigned int i = 0; i < vn.size(); i++)
            vn[i].qt[ts1] = vn[i].get_qt(sn[vn[i].i1].t[ts2] - sn[vn[i].i2].t[ts2], ts2);                               //熱量の計算
    }

    void calc_x(long ts){
        int itr = 0;
        vector<double> pre_x(sn.size(), 0.0);

        for(unsigned int i = 0; i < sn.size(); i++)                     //1ステップ前の絶対湿度を継承
            if(get<3>(sn[i].flag) == SN_CALC)   sn[i].x[ts] = sn[i].x[ts - 1];
        do{
            for(unsigned int i = 0; i < sn.size(); i++)     pre_x[i] = sn[i].x[ts]; 

            for(unsigned int i = 0; i < x_idc.size(); i++){             //節点ループ（vnループではない）
                double  sum_G_a = 0.0,                                  //自節点向けの空気流入量Ga[kg/s]を集計 
                        sum_G_m = 0.0;                                  //自節点向けの水蒸気流入量Ga[kg/s] * x[kg/kg(DA)]を集計
                for(unsigned int j = 0; j < vn.size(); j++){            //vnループ
                    if(vn[j].qv[ts] > 0 && x_idc[i] == vn[j].i2 && pre_x[vn[j].i1] != 0.0){     //自分向けで相手の絶対湿度が0.0kg/kg(DA)でない場合
                            sum_G_m += vn[j].qv[ts] * RHO20 * pre_x[vn[j].i1];
                            sum_G_a += vn[j].qv[ts] * RHO20;
                    }
                    if(vn[j].qv[ts] < 0 && x_idc[i] == vn[j].i1 && pre_x[vn[j].i2] != 0.0){
                            sum_G_m += -vn[j].qv[ts] * RHO20 * pre_x[vn[j].i2];
                            sum_G_a += -vn[j].qv[ts] * RHO20;
                    }
                }

                if(!sn[x_idc[i]].v.empty() && pre_x[x_idc[i]] != 0.0)                   //移流に伴う絶対湿度変化
                    sn[x_idc[i]].x[ts] = (sum_G_m + pre_x[x_idc[i]] * (sn[x_idc[i]].v[ts] - sum_G_a)) / sn[x_idc[i]].v[ts];
                else
                    sn[x_idc[i]].x[ts] = (sum_G_a != 0.0) ? sum_G_m / sum_G_a: 0.0;

                double mx = (sn[x_idc[i]].mx.empty()) ? 0 :sn[x_idc[i]].mx[ts];         //発湿量が設定されていれば発湿
                if(sn[x_idc[i]].v.empty())      sn[x_idc[i]].x[ts] += mx / RHO20;
                else                            sn[x_idc[i]].x[ts] += mx / (sn[x_idc[i]].v[ts] * RHO20);

            }
        }while(itr++ < sts.t_step);       //タイムステップが1hだと発散するので、1sごとに計算

    }

    int calc(){
        vector<double>  pre_p(sn.size(), 0.0), pre_t(sn.size(), 0.0);
        double          delta_p, delta_t; 

        LOG_PRINT("******************************************************************************Start calc!" << endl);     

        LOG_PRINT("length:     " << sts.length << endl);
        LOG_PRINT("sts-t_step: " << sts.t_step << endl);
        LOG_PRINT("sts-solve:  " << sts.solve << endl);
        LOG_PRINT("step_p:     " << sts.step_p << endl);
        LOG_PRINT("vent_err:   " << sts.vent_err << endl);
        LOG_PRINT("step_t:     " << sts.step_t << endl);
        LOG_PRINT("thrm_err:   " << sts.thrm_err << endl);
        LOG_PRINT("sor_ratio:  " << sts.sor_ratio << endl);
        LOG_PRINT("sor_err:    " << sts.sor_err << endl);
        LOG_CONTENTS(endl);

        for(unsigned int i = 0; i < sn.size(); i++){
            LOG_PRINT("sn[" << i << "] ="); 
            LOG_CONTENTS("(" << get<0>(sn[i].flag) << "," << get<1>(sn[i].flag)  << "," << get<2>(sn[i].flag) << "," << get<3>(sn[i].flag) << ") ");
            LOG_CONTENTS("i=" << sn[i].i << " ,s_i=" << sn[i].s_i << ", p[0]=" << sn[i].p[0] << ", c[0]" << sn[i].c[0] << ", t[0]" << sn[i].t[0] << ", x[0]" << sn[i].x[0]);
            if(sn[i].m.size()     != 0) LOG_CONTENTS(", m[0]"       << sn[i].m[0]);
            if(sn[i].h_sr.size()  != 0) LOG_CONTENTS(", h_sr[0]="   << sn[i].h_sr[0]);
            if(sn[i].h_inp.size() != 0) LOG_CONTENTS(", h_inp[0]="  << sn[i].h_inp[0]);
            if(sn[i].v.size()     != 0) LOG_CONTENTS(", v[0]="      << sn[i].v[0]);
            if(sn[i].beta.size()  != 0) LOG_CONTENTS(", beta[0]="   << sn[i].beta[0]);
            LOG_CONTENTS(endl);
        }
        for(unsigned int i = 0; i < vn.size(); i++){
            LOG_PRINT("vn[" << i << "] = " << vn[i].vn_type << " (" << vn[i].i1 << "," << vn[i].i2  << ") " << vn[i].h1[0] << " - " << vn[i].h2[0]);
            LOG_CONTENTS(", qv[0]=" << vn[i].qv[0] << ", qv[1]=" << vn[i].qv[1] << ", qt[0]=" << vn[i].qt[0] << ", qt[1]=" << vn[i].qt[1]); 
            if(vn[i].alpha.size() != 0) LOG_CONTENTS(", alpha[0]=" << vn[i].alpha[0]);
            if(vn[i].area.size()  != 0) LOG_CONTENTS(", area[0]="  << vn[i].area[0]);
            if(vn[i].a.size()     != 0) LOG_CONTENTS(", a[0]="     << vn[i].a[0]);
            if(vn[i].n.size()     != 0) LOG_CONTENTS(", n[0]="     << vn[i].n[0]);
            if(vn[i].eta.size()   != 0) LOG_CONTENTS(", eta[0]="   << vn[i].eta[0]);
            if(vn[i].q_max.size() != 0) LOG_CONTENTS(", p_max[0]=" << vn[i].q_max[0]);
            if(vn[i].p_max.size() != 0) LOG_CONTENTS(", q_max[0]=" << vn[i].p_max[0]);
            if(vn[i].q1.size()    != 0) LOG_CONTENTS(", q1[0]="    << vn[i].q1[0]);
            if(vn[i].p1.size()    != 0) LOG_CONTENTS(", p1[0]="    << vn[i].p1[0]);
            LOG_CONTENTS(endl);
        }

        for(unsigned int i = 0; i < tn.size(); i++){
            LOG_PRINT("tn[" << i << "] = " << tn[i].tn_type << " (" << tn[i].i1 << "," << tn[i].i2 << ")");
            LOG_CONTENTS(", qt[0]=" << tn[i].qt[0]);
            if(tn[i].cdtc.size()      != 0) LOG_CONTENTS(", cdtc[0]="      << tn[i].cdtc[0]);
            if(tn[i].ms.size()        != 0) LOG_CONTENTS(", ms[0]="        << tn[i].ms[0]);
            if(tn[i].area.size()      != 0) LOG_CONTENTS(", area[0]="      << tn[i].area[0]);
            if(tn[i].cof_r.size()     != 0) LOG_CONTENTS(", cof_r[0]="     << tn[i].cof_r[0]);
            if(tn[i].cof_phi.size()   != 0) LOG_CONTENTS(", cof_phi[0]="   << tn[i].cof_phi[0]);
            if(tn[i].t_dash_gs.size() != 0) LOG_CONTENTS(", t_dash_gs[0]=" << tn[i].t_dash_gs[0]);
            if(tn[i].aircon_on.size() != 0) LOG_CONTENTS(", aircon_on[0]=" << tn[i].aircon_on[0]);
            if(tn[i].ac_mode.size()   != 0) LOG_CONTENTS(", ac_mode[0]="   << tn[i].ac_mode[0]);
            if(tn[i].pre_tmp.size()   != 0) LOG_CONTENTS(", pre_tmp[0]="   << tn[i].pre_tmp[0]);
            LOG_CONTENTS(endl);
        }

        LOG_PRINT("i_vn_ac = ");
        for(unsigned int i = 0; i < i_vn_ac.size(); i++) LOG_CONTENTS(" " << i_vn_ac[i]);
        LOG_CONTENTS(endl);

        LOG_PRINT("i_tn_ac = ");
        for(unsigned int i = 0; i < i_tn_ac.size(); i++) LOG_CONTENTS(" " << i_tn_ac[i]);
        LOG_CONTENTS(endl);
        
        LOG_PRINT("t_idc.size() = " << t_idc.size() << endl);
        LOG_PRINT("c_idc.size() = " << c_idc.size() << endl);
        LOG_PRINT("v_idc.size() = " << v_idc.size() << endl); 
        LOG_PRINT("x_idc.size() = " << x_idc.size() << endl); 
        LOG_CONTENTS(endl);

        for(long ts = 1; ts < sts.length; ts++){
            LOG_PRINT("***** ts = " << ts << "******************************************************************************" << endl); 
            if(v_idc.size() > 0){
                //for(unsigned int i = 0; i < sn.size(); i++)                                                       
                //    if(get<0>(sn[i].flag) == SN_CALC)   sn[i].p[ts] = sn[i].p[ts - 1];                            //ここを実行すると、連続して換気回路網計算した際にケタ落ちエラー
                calc_qv(ts, ts - 1);
            }

            if(t_idc.size() > 0){    
                for(unsigned int i = 0; i < sn.size(); i++){
                    if(get<2>(sn[i].flag) == SN_CALC)   sn[i].t[ts] = sn[i].t[ts - 1];
                    if(get<2>(sn[i].flag) == SN_DLY)    sn[i].t[ts] = sn[sn[i].s_i].t[ts - 1];
                }
                calc_qt(ts, ts - 1);
                for(unsigned int i = 0; i < tn.size(); i++)
                    if(tn[i].tn_type == TN_GROUND)      tn[i].refresh(sn[tn[i].i1].t[ts], sn[tn[i].i2].t[ts], ts);   //地盤のリフレッシュ
            }

            do{
                delta_p = 0.0;
                delta_t = 0.0;

                for(unsigned int i = 0; i < sn.size(); i++){
                    if(get<0>(sn[i].flag) == SN_CALC)   pre_p[i] = sn[i].p[ts];
                    if(get<2>(sn[i].flag) == SN_CALC)   pre_t[i] = sn[i].t[ts];
                }
                calc_vent(ts);
                if(x_idc.size() > 0)    calc_x(ts);
                calc_thrm(ts);

                for(unsigned int i = 0; i < sn.size(); i++){
                    if(get<0>(sn[i].flag) == SN_CALC)   delta_p += pow(sn[i].p[ts] - pre_p[i], 2.0) / v_idc.size();
                    if(get<2>(sn[i].flag) == SN_CALC)   delta_t += pow(sn[i].t[ts] - pre_t[i], 2.0) / t_idc.size();
                }
                LOG_PRINT(" delta_p = " << delta_p << ", delta_t = " << delta_t << " >>>>> " << sqrt((pow(delta_p, 2.0) + pow(delta_t, 2.0))/ 2) << endl);

            }while(sts.conv_err < sqrt((pow(delta_p, 2.0) + pow(delta_t, 2.0))/ 2));
            
            if(v_idc.size() > 0){
                calc_qv(ts, ts);
                LOG_CONTENTS("p ");
                for(unsigned int i = 0; i < sn.size(); i++)    LOG_CONTENTS(sn[i].i << ": " << sn[i].p[ts] << "Pa, ");
                LOG_CONTENTS(endl << "qv ");
                for(unsigned int i = 0; i < vn.size(); i++)    LOG_CONTENTS(vn[i].i << ": " << vn[i].qv[ts] * 3600 << "m3/h, ");
                LOG_CONTENTS(endl);    
            }

            if(t_idc.size() > 0){
                calc_qt(ts, ts);
                LOG_CONTENTS("t ");
                for(unsigned int i = 0; i < sn.size(); i++)    LOG_CONTENTS(sn[i].i << ": " << sn[i].t[ts] << "deg, ");
                LOG_CONTENTS(endl << "qt1 ");
                for(unsigned int i = 0; i < vn.size(); i++)    LOG_CONTENTS(vn[i].i << ": " << vn[i].qt[ts] << "W, ");
                LOG_CONTENTS(endl << "qt2 ");
                for(unsigned int i = 0; i < tn.size(); i++)    LOG_CONTENTS(tn[i].i << ": " << tn[i].qt[ts] << "W, ");
                LOG_CONTENTS(endl);
            }

            if(x_idc.size() > 0){
                for(unsigned int i = 0; i < sn.size(); i++)    LOG_CONTENTS("x" << sn[i].i << ": " << sn[i].x[ts] << "kg/kg(DA), ");
                LOG_CONTENTS(endl);
            }

            if(c_idc.size() > 0){
                for(unsigned int i = 0; i < c_idc.size(); i++){        
                    double m;
                    if(sn[c_idc[i]].m.empty())  m = 0;
                    else                        m = sn[c_idc[i]].m[ts];

                    sn[c_idc[i]].c[ts] =  sn[c_idc[i]].c[ts - 1] * exp(-sn[c_idc[i]].beta[ts] * sts.t_step);
                    if(sn[c_idc[i]].beta[ts] == 0){
                        sn[c_idc[i]].c[ts] += m * sts.t_step / sn[c_idc[i]].v[ts];
                    }
                    else{                            
                        sn[c_idc[i]].c[ts] += m / (sn[c_idc[i]].beta[ts] * sn[c_idc[i]].v[ts]) * (1 - exp(-sn[c_idc[i]].beta[ts] * sts.t_step));
                    }
                }
                for(unsigned int i = 0; i < vn.size(); i++){
                    double eta;
                    if(vn[i].eta.empty())   eta = 0;
                    else                    eta = vn[i].eta[ts];
                    if(vn[i].qv[ts] > 0 && get<1>(sn[vn[i].i2].flag) == SN_CALC)    
                        sn[vn[i].i2].c[ts] += vn[i].qv[ts] * (sn[vn[i].i1].c[ts - 1] * (1 - eta) - sn[vn[i].i2].c[ts - 1]) * sts.t_step / sn[vn[i].i2].v[ts];
                    if(vn[i].qv[ts] < 0 && get<1>(sn[vn[i].i1].flag) == SN_CALC)    
                        sn[vn[i].i1].c[ts] += vn[i].qv[ts] * (sn[vn[i].i2].c[ts - 1] * (1 - eta) - sn[vn[i].i1].c[ts - 1]) * sts.t_step / sn[vn[i].i1].v[ts];
                }
                for(unsigned int i = 0; i < sn.size(); i++)    LOG_CONTENTS("c" << sn[i].i << ": " << sn[i].c[ts] << "num/L, ");
                LOG_CONTENTS(endl);
            }
            LOG_CONTENTS(endl);
        }
        LOG_PRINT("******************************************************************************Finish calc!" << endl << endl);
        return 0;
    }

    tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>,
          vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> result(){
        vector<vector<double>>   r_p, r_c, r_t, r_x, r_qv, r_qt1, r_qt2, r_elec;
        for(Node n: sn){
            r_p.push_back(n.p);
            r_c.push_back(n.c);
            r_t.push_back(n.t);
            r_x.push_back(n.x);
        }
        for(Vent_Net nt: vn){
            r_qv.push_back(nt.qv);
            r_qt1.push_back(nt.qt);
        }
        for(Thrm_Net nt: tn){   
            r_qt2.push_back(nt.qt);
            r_elec.push_back(nt.E_E);
        }
        return forward_as_tuple(r_p, r_c, r_t, r_x, r_qv, r_qt1, r_qt2, r_elec);
    }
};
