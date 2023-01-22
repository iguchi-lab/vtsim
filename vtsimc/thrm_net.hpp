#include "RAC_model.hpp"
#include "Duct_Central_model.hpp"

using namespace std;

class Thrm_Net{
public:
    int i, i1, i2, tn_type;                                                                 //入口出口のノード番号、熱回路網の種類
    
    vector<double>      cdtc;                                                               //コンダクタンス
    vector<double>      ms;                                                                 //日射取得率
    vector<double>      area;                                                               //面積
    double              phi_0;                                                              //応答係数                                                        
    vector<double>      cof_r, cof_phi, t_dash_gs;
    vector<double>      qt;                                                                 //熱流
    vector<double>      Ls = {0.0}, Ll = {0.0}, E_E = {0.0};                                //顕熱負荷、潜熱負荷、電力

    int                 ac_model;
    vector<int>         aircon_on, ac_mode;                                                 //エアコンのON/OFF、エアコン運転モード
    vector<double>      pre_tmp, pre_x, To, ho;                                             //エアコン設定温度
    double              q_rtd_C, q_max_C, e_rtd_C, q_rtd_H, q_max_H, e_rtd_H;
    double              q_hs_rtd_H, q_hs_mid_H, q_hs_min_H, P_hs_rtd_H, 
                        V_fan_rtd_H, V_fan_mid_H, P_fan_rtd_H, P_fan_mid_H, V_hs_dsgn_H,
                        q_hs_rtd_C, q_hs_mid_C, q_hs_min_C, P_hs_rtd_C, 
                        V_fan_rtd_C, V_fan_mid_C, P_fan_rtd_C, P_fan_mid_C, V_hs_dsgn_C, 
                        Theta, h;
    bool                dualcompressor;             

    Thrm_Net(long length, int i_, int i1_ , int i2_, int tn_type_){
        i       = i_;
        i1      = i1_;
        i2      = i2_;
        tn_type = tn_type_;
        qt.assign(length, 0.0);
        t_dash_gs.assign(10, 0.0);
        aircon_on.assign(length, 0);                                                        //熱量を0に初期化      
        E_E.assign(length, 0.0);    
    }              

    double get_qt(double dt, long ts){
        double sum_t_dash_gs = 0.0;                                                         //項別成分の合計
        switch(tn_type){
            case TN_SIMPLE:    
                return cdtc[ts] * dt;                                                       //単純熱回路網の場合
            case TN_GROUND:                                                                 //地盤の場合  
                for(int i = 0; i < 10; i++)     sum_t_dash_gs += t_dash_gs[i];              //合計の計算
                return area[ts] / phi_0 * (dt - sum_t_dash_gs);                             //熱流の計算
            default:
                return 0.0;
        }
    }

    double get_qt_s(double h_sr, long ts)  {return ms[ts] * h_sr;}                          //日射取得
    double get_qt_h(double h_inp, long ts) {return h_inp;}                                  //発熱
    
    void refresh(double t_uf_t_1, double t_g_t_1, long ts){
        double sum_t_dash_gs = 0.0;                                                         //項別成分の合計
        for(int i = 0; i < 10; i++)         sum_t_dash_gs += t_dash_gs[i];                  //合計の計算
        for(int i = 0; i < 10; i++)
            t_dash_gs[i] = cof_phi[i] * qt[ts] / area[ts] + cof_r[i] * t_dash_gs[i];        //吸熱応答の項別成分の計算
    }

    /*********************************************************************************************
     電力_RAC
    **********************************************************************************************/
    tuple<double , double> get_RAC_hs_out(double Theta_in, double Theta_out, double X_in, double X_out,
                                      double V_supply, double Theta, double h, double ts){

        double  L_H   = c_p_air * rho_air * (Theta_out - Theta_in ) * V_supply * 3600 / 1e+6,
                L_CS  = c_p_air * rho_air * (Theta_in  - Theta_out) * V_supply * 3600 / 1e+6, 
                L_CL  = L_wtr   * rho_air * (X_in      - X_out) * V_supply * 1e+3 * 3600 / 1e+6;

        if((L_H > 0.0) && (L_CS < 0.0)){                          //暖房時
            double  q_r_max_H, Q_r_max_H, Q_max_H, Q_T_H, Q_dash_T_H;

            if((Theta < 5.0) && (h >= 80.0))    C_df_H = 0.77;                      //外気温5.0℃未満、相対湿度80%以上はデフロスト C_df_H=0.77
            else                                C_df_H = 1.00;

            q_r_max_H  = q_max_H / q_rtd_H;                                         //(4) 定格暖房能力に対する最大暖房能力の比 [-]
            Q_r_max_H  = calc_Q_r_max_H(q_rtd_C, q_r_max_H, Theta);                 //定格暖房能力に対する最大暖房出力の比 [-]
            Q_max_H    = Q_r_max_H * q_rtd_H * C_af_H * C_df_H * 3600 * 1e-6;       //(1) １時間当たりの最大暖房出力 [MJ/h] 
            Q_T_H      = min(Q_max_H, L_H);                                         //１時間当たりの暖房設備機器の処理暖房負荷 [MJ/h]
            Q_dash_T_H = Q_T_H * (1.0 / (C_af_H * C_df_H));                         //補正処理暖房負荷 [MJ/h]

            Theta_out = Theta_in + Q_T_H / (c_p_air * rho_air * V_supply * 3600 / 1e+6);
            
            if(Q_dash_T_H > 0.0){
                double x1, x2;
                x1 = calc_f_H_Theta(Q_dash_T_H / (q_max_H * 3600 * 1e-6), Theta, q_rtd_C, dualcompressor);
                x2 = calc_f_H_Theta(1.0 / q_r_max_H,                      7.0,   q_rtd_C, dualcompressor);
                E_E[ts] =  x1 / x2 * (q_rtd_H / e_rtd_H) * 1e-3; //(5)
                Ls[ts]  = c_p_air * rho_air * (Theta_out - Theta_in ) * V_supply        * 3600 / 1e+6;
                Ll[ts]  = L_wtr   * rho_air * (X_out     - X_in)      * V_supply * 1e+3 * 3600 / 1e+6;
            }
            return make_tuple(Theta_out, X_in);
        }        
        else if((L_H < 0.0) && (L_CS > 0.0)){
            double  q_r_max_C, Q_r_max_C, Q_max_C, L_max_CL, L_dash_CL, L_dash_C, 
                    SHF_dash, Q_max_CS, Q_max_CL, Q_T_CS, Q_T_CL, Q_dash_T_C;

            q_r_max_C  = q_max_C / q_rtd_C;                                         //(14) 定格冷房能力に対する最大冷房能力の比 [-]
            Q_r_max_C  = calc_Q_r_max_C(q_r_max_C, q_rtd_C, Theta);                 //定格冷房能力に対する最大冷房出力の比 [-]
            Q_max_C    = Q_r_max_C * q_rtd_C * C_af_C * C_hm_C * 3600 * 1e-6;       //(11) １時間当たりの最大冷房出力 [MJ/h]
            L_max_CL   = L_CS * ((1.0 - SHF_L_min_c) / SHF_L_min_c);                //(19) 最大で処理できる冷房潜熱負荷 [MJ/h]
            L_dash_CL  = min(L_max_CL, L_CL);                                       //(18) １時間当たりの補正冷房潜熱負荷 [MJ/h]
            L_dash_C   = L_CS + L_dash_CL;                                          //(17) １時間当たりの補正冷房負荷 [MJ/h]
            SHF_dash   = (L_dash_C != 0) ? L_CS / L_dash_C: 0.0;                    //(6) 冷房負荷補正顕熱比 [-]       
            Q_max_CS   = Q_max_C * SHF_dash;                                        //(15a) １時間当たりの最大冷房顕熱出力 [MJ/h]
            Q_max_CL   = min(Q_max_C * (1.0 - SHF_dash), L_dash_CL);                //(15b) １時間当たりの最大冷房潜熱出力 [MJ/h]
            Q_T_CS     = min(Q_max_CS, L_CS);                                       //1時間当たりの冷房設備機器の処理冷房顕熱負荷 [MJ/h] 
            Q_T_CL     = min(Q_max_CL, L_CL);                                       //1時間当たりの冷房設備機器の処理冷房潜熱負荷 [MJ/h] 
            Q_dash_T_C = (Q_T_CS + Q_T_CL) * (1.0 / (C_hm_C * C_af_C));             //(25) 1時間当たりの冷房設備機器の処理冷房顕熱負荷 [MJ/h] 

            Theta_out = Theta_in - Q_T_CS / (c_p_air * rho_air * V_supply * 3600 / 1e+6);
            X_out     = X_in     - Q_T_CL / (L_wtr   * rho_air * V_supply * 1e+3 * 3600 / 1e+6);

            if(Q_dash_T_C > 0.0){
                double x1, x2;
                x1 = calc_f_C_Theta(Q_dash_T_C / (q_max_C * 3600 * 1e-6), Theta, q_rtd_C, dualcompressor);
                x2 = calc_f_C_Theta(1.0 / q_r_max_C,                        35.0, q_rtd_C, dualcompressor);
                E_E[ts] = x1 / x2 * (q_rtd_C / e_rtd_C) * 1e-3; //(20)
                Ls[ts]  = c_p_air * rho_air * (Theta_out - Theta_in ) * V_supply        * 3600 / 1e+6;
                Ll[ts]  = L_wtr   * rho_air * (X_out     - X_in)      * V_supply * 1e+3 * 3600 / 1e+6;
            }
        }
        return make_tuple(Theta_out, X_out);
    }

    /*********************************************************************************************
     電力_DUCT_CENTRAL
    **********************************************************************************************/
    tuple<double , double> get_DUCT_C_hs_out(double Theta_in, double Theta_out, double X_in, double X_out,
                                             double V_supply, double Theta, double h, double ts){
    
        double E_E_comp, E_E_fan;

        if(Theta_out > Theta_in){                          //暖房時
            tie(Theta_out, X_out, E_E_comp, E_E_fan) = calc_DUCT_E_E_H(q_hs_rtd_H, q_hs_mid_H, q_hs_min_H, P_hs_rtd_H, V_fan_rtd_H * 3600, 
                                                                       V_fan_mid_H * 3600, P_fan_rtd_H, P_fan_mid_H, V_hs_dsgn_H * 3600,
                                                                       Theta, h, Theta_out, Theta_in, X_out, X_in, V_supply * 3600);
            E_E[ts] = E_E_comp + E_E_fan;
            Ls[ts]  = c_p_air * rho_air * (Theta_out - Theta_in ) * V_supply        * 3600 / 1e+6;
            Ll[ts]  = L_wtr   * rho_air * (X_out     - X_in)      * V_supply * 1e+3 * 3600 / 1e+6;
            return make_tuple(Theta_out, X_out);
        }
        else if(Theta_out < Theta_in){                      //冷房時
            tie(Theta_out, X_out, E_E_comp, E_E_fan) = calc_DUCT_E_E_C(q_hs_rtd_C, q_hs_mid_C, q_hs_min_C, P_hs_rtd_C, V_fan_rtd_C * 3600, 
                                                                       V_fan_mid_C * 3600, P_fan_rtd_C, P_fan_mid_C, V_hs_dsgn_C * 3600,
                                                                       Theta, h, Theta_out, Theta_in, X_out, X_in, V_supply * 3600);
            E_E[ts] = E_E_comp + E_E_fan;
            Ls[ts]  = c_p_air * rho_air * (Theta_out - Theta_in ) * V_supply        * 3600 / 1e+6;
            Ll[ts]  = L_wtr   * rho_air * (X_out     - X_in)      * V_supply * 1e+3 * 3600 / 1e+6;
            return make_tuple(Theta_out, X_out);                                         
        }
        return make_tuple(0.0, 0.0);
    }
};