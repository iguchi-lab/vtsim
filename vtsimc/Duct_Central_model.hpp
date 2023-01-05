#include <iostream>
#include <tuple>
#include <math.h>

using namespace::std;

//飽和蒸気の温度から圧力を求める関数 (13)
double get_f_p_sgas(double Theta){
    return   2.75857926950901 * 1e-17 * pow(Theta, 8) \
           + 1.49382057911753 * 1e-15 * pow(Theta, 7) \
           + 6.52001687267015 * 1e-14 * pow(Theta, 6) \
           + 9.14153034999975 * 1e-12 * pow(Theta, 5) \
           + 3.18314616500361 * 1e-9  * pow(Theta, 4) \
           + 1.60703566663019 * 1e-6  * pow(Theta, 3) \
           + 3.06278984019513 * 1e-4  * pow(Theta, 2) \
           + 2.54461992992037 * 1e-2  * Theta \
           + 7.98086455154775 * 1e-1;
}

//圧縮機吸引領域において過熱蒸気の圧力と温度から比エンタルピーを求める関数 (14)
double get_f_H_gas_comp_in(double P, double Theta){
    double K = Theta + 273.15;
    double K2 = K  * K;
    double K3 = K2 * K;
    double P2 = P  * P;
    double P3 = P2 * P;
    //double P4 = P2 * P2;
    return -1.00110355   * 1e-1 * P3 \
           - 1.184450639 * 10   * P2 \
           - 2.052740252 * 1e+2 * P \
           + 3.20391     * 1e-6 * K3 \
           - 2.24685     * 1e-3 * K2 \
           + 1.279436909        * K \
           + 3.1271238   * 1e-2 * P2 * K \
           - 1.415359    * 1e-3 * P  * K2 \
           + 1.05553912         * P  * K \
           + 1.949505039 * 1e+2;
}

//圧縮機吐出領域において過熱蒸気の圧力と比エントロピーから比エンタルピーを求める関数 (15)
double get_f_H_gas_comp_out(double P, double S){
    double P2 = P  * P;
    double P3 = P2 * P;
    double P4 = P2 * P2;
    double S2 = S  * S;
    double S3 = S2 * S;
    double S4 = S2 * S2;
    return  -1.869892835947070 * 1e-1 * P4 \
           + 8.223224182177200 * 1e-1 * P3 \
           + 4.124595239531860        * P2 \
           - 8.346302788803210 * 1e+1 * P \
           - 1.016388214044490 * 1e+2 * S4 \
           + 8.652428629143880 * 1e+2 * S3 \
           - 2.574830800631310 * 1e+3 * S2 \
           + 3.462049327009730 * 1e+3 * S \
           + 9.209837906396910 * 1e-1 * P3 * S \
           - 5.163305566700450 * 1e-1 * P2 * S2 \
           + 4.076727767130210 * P         * S3 \
           - 8.967168786520070 * P2        * S \
           - 2.062021416757910 * 10   * P  * S2 \
           + 9.510257675728610 * 10   * P  * S \
           - 1.476914346214130 * 1e+3;
}

//過熱蒸気の圧力と比エンタルピーから比エントロピーを求める関数 (16)
double get_f_S_gas(double P, double h){
    double P2 = P  * P;
    double P3 = P2 * P;
    double P4 = P2 * P2;
    double h2 = h  * h;
    double h3 = h2 * h;
    double h4 = h2 * h2;

    return   5.823109493752840 * 1e-2 * P4 \
           - 3.309666523931270 * 1e-1 * P3 \
           + 7.700179914440890 * 1e-1 * P2 \
           - 1.311726004718660        * P \
           + 1.521486605815750 * 1e-9 * h4 \
           - 2.703698863404160 * 1e-6 * h3 \
           + 1.793443775071770 * 1e-3 * h2 \
           - 5.227303746767450 * 1e-1 * h \
           + 1.100368875131490 * 1e-4 * pow(P, 3) * h \
           + 5.076769807083600 * 1e-7 * pow(P, 2) * h2 \
           + 1.202580329499520 * 1e-8 * P * h3 \
           - 7.278049214744230 * 1e-4 * pow(P, 2) * h \
           - 1.449198550965620 * 1e-5 * P * h2 \
           + 5.716086851760640 * 1e-3 * P * h \
           + 5.818448621582900 * 1e+1;
}

//過冷却液の圧力と温度から比エンタルピーを求める関数 (17)
double get_f_H_liq(double P, double Theta){
    double K  = Theta + 273.15;
    double K2 = K  * K;
    double K3 = K2 * K;
    double P2 = P  * P;
    double P3 = P2 * P;

    return   1.7902915   * 1e-2 * P3 \
           + 7.96830322  * 1e-1 * P2 \
           + 5.985874958 * 10   * P \
           + 0                  * K3 \
           + 9.86677     * 1e-4 * K2 \
           + 9.8051677   * 1e-1 * K \
           - 3.58645     * 1e-3 * pow(P, 2) * K \
           + 8.23122     * 1e-4 * P * K2 \
           - 4.42639115  * 1e-1 * P * K \
           - 1.415490404 * 1e+2;
}

double calc_e_ref_H_th(double Theta_ref_evp, double Theta_ref_cnd, double Theta_ref_SC, double Theta_ref_SH){
    double P_ref_evp, P_ref_cnd, Theta_ref_cnd_out, h_ref_cnd_out, 
           Theta_ref_comp_in, P_ref_comp_in, h_ref_comp_in, S_ref_comp_in, 
           S_ref_comp_out, P_ref_comp_out, h_ref_comp_out, e_ref_H_th;

    P_ref_evp         = get_f_p_sgas(Theta_ref_evp);                                            //蒸発圧力 (12)
    P_ref_cnd         = get_f_p_sgas(Theta_ref_cnd);                                            //凝縮圧力 (11)
    Theta_ref_cnd_out = Theta_ref_cnd - Theta_ref_SC;                                           //凝縮器出力温度 (10)
    h_ref_cnd_out     = get_f_H_liq(P_ref_cnd, Theta_ref_cnd_out);                              //凝縮器出口比エンタルピー (9)
    Theta_ref_comp_in = Theta_ref_evp + Theta_ref_SH;                                           //圧縮機吸込温度 (8)
    P_ref_comp_in     = P_ref_evp;                                                              //圧縮機吸込圧力 (7)
    h_ref_comp_in     = get_f_H_gas_comp_in(P_ref_comp_in, Theta_ref_comp_in);                  //圧縮機吸込エンタルピー (6)
    S_ref_comp_in     = get_f_S_gas(P_ref_comp_in, h_ref_comp_in);                              //圧縮機吸込比エントロピー (5)
    S_ref_comp_out    = S_ref_comp_in;                                                          //圧縮機吐出比エントロピー (4)  
    P_ref_comp_out    = P_ref_cnd;                                                              //圧縮機吐出圧力 (3)       
    h_ref_comp_out    = get_f_H_gas_comp_out(P_ref_comp_out, S_ref_comp_out);                   //圧縮機吐出比エンタルピー (2)
    e_ref_H_th        = (h_ref_comp_out - h_ref_cnd_out) / (h_ref_comp_out - h_ref_comp_in);    //ヒートポンプサイクルの理論暖房効率 (1)
    return e_ref_H_th;
}

tuple<double, double, double, double> calc_DUCT_E_E_H(double q_hs_rtd_H, double q_hs_mid_H, double q_hs_min_H, double P_hs_rtd_H, 
                                                      double V_fan_rtd_H, double V_fan_mid_H, double P_fan_rtd_H, double P_fan_mid_H, double V_hs_dsgn_H,
                                                      double Theta, double h, double Theta_hs_out, double Theta_hs_in, double X_hs_out, double X_hs_in, 
                                                      double V_hs_supply){
    
    double  C_df_H = 1.0, q_hs_H, alpha_c_hex_H,
            Theta_sur_f_hex_H, Theta_ref_cnd_H, Theta_ref_evp_H, Theta_ref_SC_H, Theta_ref_SH_H,
            e_dash_th_mid_H, e_th_mid_H, e_th_rtd_H, e_th_H, e_hs_rtd_H, e_r_rtd_H, e_r_min_H, e_r_mid_H, e_r_H, e_hs_H,
            E_E_comp_H = 0.0, E_E_fan_H = 0.0;
    
    C_df_H    = ((Theta < 5.0) && (h >= 80.0)) ? 0.77: 1.0;
    q_hs_H    = c_p_air * rho_air * (Theta_hs_out - Theta_hs_in) * (V_hs_supply / 3600) * (1 / C_df_H);                                 //(3) 1時間当たりの熱源機の平均暖房能力 [W]
    E_E_fan_H = (q_hs_H > 0) ? P_fan_rtd_H * V_hs_supply / V_hs_dsgn_H * 1e-3 : 0.0;                                                    //(37) 1時間当たりの送風機の消費電力量のうちの暖房設備への付加分 [kWh/h] 

    /*中間*************************************************************************************************************************/
    alpha_c_hex_H     = (-0.0017 * pow((V_fan_mid_H / 3600) / A_f_hex, 2) + 0.044 * ((V_fan_mid_H / 3600) / A_f_hex) + 0.0271) * 1e+3;  //(35) 暖房時の室内熱交換器表面の顕熱伝達率 [-]
    Theta_sur_f_hex_H = 20 + (q_hs_mid_H / (2 * c_p_air * rho_air * V_fan_mid_H)) * 3600 + (q_hs_mid_H / (A_e_hex * alpha_c_hex_H));    //(33) 暖房時の室内機熱交換器の表面温度 [℃]
    Theta_ref_cnd_H   = (Theta_sur_f_hex_H > 65) ? 65: Theta_sur_f_hex_H;                                                               //(23) 暖房時の冷媒の凝縮温度 [℃]        
    Theta_ref_evp_H   = min(max(7.0 - (0.100 * Theta_ref_cnd_H + 2.95), -50.0), Theta_ref_cnd_H - 5.0);                                 //(24) 暖房時の冷媒の蒸発温度 [℃]
    Theta_ref_SC_H    = 0.245 * Theta_ref_cnd_H - 1.72;                                                                                 //(25) 暖房時の冷媒の過冷却度 [℃]
    Theta_ref_SH_H    = 4.49 - 0.036 * Theta_ref_cnd_H;                                                                                 //(26) 暖房時の冷媒の過熱度 [℃]
    e_dash_th_mid_H   = calc_e_ref_H_th(Theta_ref_evp_H, Theta_ref_cnd_H, Theta_ref_SC_H, Theta_ref_SH_H);                              //4_8_a (1) ヒートポンプサイクルの理論暖房効率
    e_th_mid_H        = e_dash_th_mid_H;                                                                                                //(20) 中間暖房能力運転時のヒートポンプサイクルの理論効率 [-];                                                                                                
   /*定格**************************************************************************************************************************/
    alpha_c_hex_H     = (-0.0017 * pow((V_fan_rtd_H / 3600) / A_f_hex,  2) + 0.044 * ((V_fan_rtd_H / 3600) / A_f_hex) + 0.0271) * 1e+3; //(35) 暖房時の室内熱交換器表面の顕熱伝達率 [-]
    Theta_sur_f_hex_H = 20 + (q_hs_rtd_H / (2 * c_p_air * rho_air * V_fan_rtd_H)) * 3600 + (q_hs_rtd_H / (A_e_hex * alpha_c_hex_H));    //(33) 暖房時の室内機熱交換器の表面温度 [℃]       
    Theta_ref_cnd_H   = (Theta_sur_f_hex_H > 65) ? 65: Theta_sur_f_hex_H;                                                               //(23) 暖房時の冷媒の凝縮温度 [℃]        
    Theta_ref_evp_H   = min(max(7.0 - (0.100 * Theta_ref_cnd_H + 2.95), -50.0), Theta_ref_cnd_H - 5.0);                                 //(24) 暖房時の冷媒の蒸発温度 [℃]
    Theta_ref_SC_H    = 0.245 * Theta_ref_cnd_H - 1.72;                                                                                 //(25) 暖房時の冷媒の過冷却度 [℃]
    Theta_ref_SH_H    = 4.49 - 0.036 * Theta_ref_cnd_H;                                                                                 //(26) 暖房時の冷媒の過熱度 [℃]
    e_th_rtd_H        = calc_e_ref_H_th(Theta_ref_evp_H, Theta_ref_cnd_H, Theta_ref_SC_H, Theta_ref_SH_H);                              //(19) 定格暖房能力運転時のヒートポンプサイクルサイクルの理論効率 [-]
   /******************************************************************************************************************************/
    alpha_c_hex_H     = (-0.0017 * pow((V_hs_supply / 3600) / A_f_hex, 2) + 0.044 * ((V_hs_supply / 3600) / A_f_hex) + 0.0271) * 1e+3;  //(35) 暖房時の室内熱交換器表面の顕熱伝達率 [-]
    Theta_sur_f_hex_H = ((Theta_hs_in + Theta_hs_out) / 2) + (c_p_air * rho_air * V_hs_supply * (Theta_hs_out - Theta_hs_in) / 3600) 
                        / (A_e_hex * alpha_c_hex_H);                                                                                    //(31) 暖房時の室内機熱交換器の表面温度 [℃] 
    Theta_ref_cnd_H   = (Theta_sur_f_hex_H > 65) ? 65: Theta_sur_f_hex_H;                                                               //(23) 暖房時の冷媒の凝縮温度 [℃]
    Theta_ref_evp_H   = min(max(Theta - (0.100 * Theta_ref_cnd_H + 2.95), -50.0), Theta_ref_cnd_H - 5.0);                               //(24) 暖房時の冷媒の蒸発温度 [℃]
    Theta_ref_SC_H    = 0.245 * Theta_ref_cnd_H - 1.72;                                                                                 //(25) 暖房時の冷媒の過冷却度 [℃]
    Theta_ref_SH_H    = 4.49 - 0.036 * Theta_ref_cnd_H;                                                                                 //(26) 暖房時の冷媒の過熱度 [℃]
    e_th_H            = calc_e_ref_H_th(Theta_ref_evp_H, Theta_ref_cnd_H, Theta_ref_SC_H, Theta_ref_SH_H);                              //(17) 日付dの時刻tにおける暖房時のヒートポンプサイクルの理論効率 [-]
    /*****************************************************************************************************************************/
    e_hs_rtd_H        = q_hs_rtd_H / (P_hs_rtd_H - P_fan_rtd_H);                                                                        //(11) 定格暖房能力運転時のヒートポンプサイクルの理論効率に対する熱源機の効率の比 [-]
    e_r_rtd_H         = min(max(e_hs_rtd_H / e_th_rtd_H, 0.0), 1.0);    
    e_r_min_H         = e_r_rtd_H * 0.65;                                                                                               //(15) 最小暖房能力運転時のヒートポンプサイクルの理論効率に対する熱源機の効率の比 [-]
    e_r_mid_H         = min(max(e_r_rtd_H * 0.95,  0.0), 1.0);                                                                          //(13) 中間暖房能力運転時のヒートポンプサイクルの理論効率に対する熱源機の効率の比
    //e_hs_mid_H = q_hs_mid_H / (P_hs_mid_H - P_fan_mid_H)    #定格能力試験と中間能力試験の値を入力する場合はこちら
    //e_r_mid_H = e_hs_mid_H / e_th_mid_H

    //(9) 暖房時のヒートポンプサイクルの理論効率に対する熱源機の効率の比 [-]
    if      (q_hs_H      <= q_hs_min_H)                         e_r_H = e_r_min_H - (q_hs_min_H - q_hs_H) * (e_r_min_H / q_hs_min_H);
    else if ((q_hs_min_H < q_hs_H) && (q_hs_H <= q_hs_mid_H))   e_r_H = e_r_mid_H - (q_hs_mid_H - q_hs_H) * ((e_r_mid_H - e_r_min_H) / (q_hs_mid_H - q_hs_min_H));
    else if ((q_hs_mid_H < q_hs_H) && (q_hs_H <= q_hs_rtd_H))   e_r_H = e_r_rtd_H - (q_hs_rtd_H - q_hs_H) * ((e_r_rtd_H - e_r_mid_H) / (q_hs_rtd_H - q_hs_mid_H));
    else if ((q_hs_rtd_H < q_hs_H) && (e_r_rtd_H > 0.4))        e_r_H = max(e_r_rtd_H - (q_hs_H - q_hs_rtd_H) * (e_r_rtd_H / q_hs_rtd_H), 0.4);
    else if ((q_hs_rtd_H < q_hs_H) && (e_r_rtd_H <= 0.4))       e_r_H = e_r_rtd_H;
    else                                                        e_r_H = 0.0;

    e_hs_H = e_th_H * e_r_H;                                                                                                            //(7) 日付dの時刻tにおける暖房時の熱源機の効率 [-]                                                                                                    
    E_E_comp_H = ((q_hs_H > 0) && (e_hs_H > 0.0)) ? (q_hs_H / e_hs_H) * 1e-3 : 0.0;                                                                         //(5)
    
    double Theta_out = Theta_hs_in + q_hs_H / (c_p_air * rho_air * V_hs_supply * 3600 / 1e+6);
    
    return make_tuple(Theta_out, X_hs_in, E_E_comp_H, E_E_fan_H);
}

//(34) 冷房時の室内機熱交換器の表面温度 [℃]
double  func(double x, double V_fan_x_C, double alpha_dash_c_hex_C, double alpha_c_hex_C, double q_hs_x_C){
    double  q_hs_CS, T, k, P_vs, X_surf_hex_C, q_hs_CL,
            a1 = -6096.9385,    a2 = 21.2409642,    a3 = -0.02711193,   a4 = 0.00001673952,     a5 = 2.433502,
            b1 = -6024.5282,    b2 = 29.32707,      b3 = 0.010613863,   b4 = -0.000013198825,   b5 = -0.49382577;

    q_hs_CS = (27 - x) / (3600 / (2 * c_p_air * rho_air * V_fan_x_C) + 1 / (A_e_hex * alpha_c_hex_C));
    T       = x + 273.16;                                                                                                               //(6) 絶対温度 [K]
    k = (x > 0) ? a1 / T + a2 + a3 * T + a4 * pow(T, 2) + a5 * log(T): b1 / T + b2 + b3 * T + b4 * pow(T, 2) + b5 * log(T);             //(5b)
    P_vs = exp(k);                                                                                                                      //(5a) 飽和水蒸気圧 [Pa]
    X_surf_hex_C = 0.622 * (P_vs / (F - P_vs));                                                                                         //(3) 飽和空気の絶対湿度 [kg/kg']
    q_hs_CL = max((0.010376 - X_surf_hex_C) / (3600 / (2 * L_wtr * rho_air * V_fan_x_C * 1e+3) + 1 / (L_wtr * A_e_hex * alpha_dash_c_hex_C * 1e+3)), 0.0);

    return q_hs_CS + q_hs_CL - q_hs_x_C;
}

tuple<double, double, double, double > calc_DUCT_E_E_C(double q_hs_rtd_C, double q_hs_mid_C, double q_hs_min_C, double P_hs_rtd_C, 
                                                       double V_fan_rtd_C, double V_fan_mid_C, double P_fan_rtd_C, double P_fan_mid_C, double V_hs_dsgn_C,
                                                       double Theta, double h, double Theta_hs_out, double Theta_hs_in, double X_hs_out, double X_hs_in, 
                                                       double V_hs_supply){
    
    double  q_hs_CS, q_hs_CL, q_hs_C, alpha_dash_c_hex_C, alpha_c_hex_C,
            Theta_sur_f_hex_C, Theta_ref_evp_C, Theta_ref_cnd_C, Theta_ref_SC_C, Theta_ref_SH_C,
            e_dash_th_mid_C, e_th_mid_C, e_dash_th_rtd_C, e_th_rtd_C, e_dash_th_C, e_th_C, e_hs_rtd_C, e_r_rtd_C, e_r_min_C, e_r_mid_C, e_r_C, e_hs_C,
            E_E_comp_C = 0.0, E_E_fan_C = 0.0;

    double tmp_t = 0.0, step_t = 1e-6, err0 = 999.0, err1;
    int i = 0;

    //cout << q_hs_rtd_C << ", " << q_hs_mid_C << ", " << q_hs_min_C << ", " << P_hs_rtd_C << ", " << 
    //       V_fan_rtd_C << ", " << V_fan_mid_C << ", " << P_fan_rtd_C << ", " << P_fan_mid_C << ", " << V_hs_dsgn_C << ", " << 
    //        Theta << ", " <<  h << ", " << Theta_hs_out << ", " << Theta_hs_in << ", " << X_hs_out << ", " << X_hs_in << ", " << V_hs_supply << endl;

    q_hs_CS = c_p_air * rho_air * (Theta_hs_in - Theta_hs_out) * (V_hs_supply / 3600);                                                  //(4b-2) 1時間当たりの熱源機の平均冷房顕熱能力 [W] 
    q_hs_CL = L_wtr * rho_air * (X_hs_in - X_hs_out) * (V_hs_supply / 3600) * 1e+3;                                                     //(4c-2) 1時間当たりの熱源機の平均冷房潜熱能力 [W]
    q_hs_C  = q_hs_CS + q_hs_CL;                                                                                                        //(4a-2) 1時間当たりの熱源機の平均冷房能力 [W]
    E_E_fan_C = (q_hs_C > 0) ? P_fan_rtd_C * V_hs_supply / V_hs_dsgn_C * 1e-3 : 0.0;                                                    //(38) 1時間当たりの送風機の消費電力量 [kWh/h] 
    //cout << q_hs_CS << ", " << q_hs_CL << ", " << q_hs_C << ", " << E_E_fan_C << endl;    
    /*****************************************************************************************************************************/    
    alpha_dash_c_hex_C = 0.050 * log((max(V_fan_mid_C, 400.0) / 3600) / A_f_hex) + 0.073;                                               //(36b) 冷房時の室内熱交換器表面の潜熱伝達率 [kg/(m2・s)] 
    alpha_c_hex_C = alpha_dash_c_hex_C * (c_p_air + c_p_w * 0.010376);                                                                  //(36a) 冷房時の室内熱交換器表面の顕熱伝達率 [W/(m2・K)]

    tmp_t = 0.0;
    err0  = 999.0;
    i = 0;
    while(abs(err0) > 1.0e-9){
        err0 = func(tmp_t,          V_fan_mid_C, alpha_dash_c_hex_C, alpha_c_hex_C, q_hs_mid_C);
        err1 = func(tmp_t + step_t, V_fan_mid_C, alpha_dash_c_hex_C, alpha_c_hex_C, q_hs_mid_C);
        tmp_t += err0 / ((err0 - err1) / step_t);
        if(i > 100){
            return make_tuple(0.0, 0.0, 0.0, 0.0);
        }
    }
    Theta_sur_f_hex_C = tmp_t;
    Theta_ref_evp_C = max(Theta_sur_f_hex_C, -50.0);                                                                                    //(28) 冷房時の冷媒の蒸発温度 [℃]
    Theta_ref_cnd_C = max(min(max(35 + 27.4 - 1.35 * Theta_ref_evp_C, 35.0), 65.0), Theta_ref_evp_C + 5.0);                             //(27) 冷房時の冷媒の凝縮温度 [℃]
    Theta_ref_SC_C = max(0.772 * Theta_ref_cnd_C - 25.6, 0.0);                                                                          //(29) 冷房時の冷媒の過冷却度 [℃]
    Theta_ref_SH_C = max(0.194 * Theta_ref_cnd_C - 3.86, 0.0);                                                                          //(30) 冷房時の冷媒の過熱度 [℃]
    e_dash_th_mid_C = calc_e_ref_H_th(Theta_ref_evp_C, Theta_ref_cnd_C, Theta_ref_SC_C, Theta_ref_SH_C);                                //# 4_8_a (1) ヒートポンプサイクルの理論暖房効率 [-]
    e_th_mid_C = e_dash_th_mid_C - 1;                                                                                                   //(22) 中間暖房能力運転時のヒートポンプサイクルの理論効率 [-]
    /*****************************************************************************************************************************/ 
    alpha_dash_c_hex_C = 0.050 * log((max(V_fan_rtd_C, 400.0) / 3600) / A_f_hex) + 0.073;                                               //(36b) 冷房時の室内熱交換器表面の潜熱伝達率 [kg/(m2・s)] 
    alpha_c_hex_C = alpha_dash_c_hex_C * (c_p_air + c_p_w * 0.010376);                                                                  //(36a) 冷房時の室内熱交換器表面の顕熱伝達率 [W/(m2・K)]

    tmp_t = 0.0;
    err0  = 999.0;
    i = 0;
    while(abs(err0) > 1.0e-9){
        err0 = func(tmp_t,          V_fan_rtd_C, alpha_dash_c_hex_C, alpha_c_hex_C, q_hs_rtd_C);
        err1 = func(tmp_t + step_t, V_fan_rtd_C, alpha_dash_c_hex_C, alpha_c_hex_C, q_hs_rtd_C);
        tmp_t += err0 / ((err0 - err1) / step_t);
        if(i > 100){
            return make_tuple(0.0, 0.0, 0.0, 0.0);
        }
    }
    Theta_sur_f_hex_C = tmp_t;
    Theta_ref_evp_C = max(Theta_sur_f_hex_C, -50.0);                                                                                    //(28) 冷房時の冷媒の蒸発温度 [℃]
    Theta_ref_cnd_C = max(min(max(35 + 27.4 - 1.35 * Theta_ref_evp_C, 35.0), 65.0), Theta_ref_evp_C + 5.0);                             //(27) 冷房時の冷媒の凝縮温度 [℃]
    Theta_ref_SC_C = max(0.772 * Theta_ref_cnd_C - 25.6, 0.0);                                                                          //(29) 冷房時の冷媒の過冷却度 [℃]
    Theta_ref_SH_C = max(0.194 * Theta_ref_cnd_C - 3.86, 0.0);                                                                          //(30) 冷房時の冷媒の過熱度 [℃]
    e_dash_th_rtd_C = calc_e_ref_H_th(Theta_ref_evp_C, Theta_ref_cnd_C, Theta_ref_SC_C, Theta_ref_SH_C);                                //# 4_8_a (1) ヒートポンプサイクルの理論暖房効率 [-]
    e_th_rtd_C = e_dash_th_rtd_C - 1;                                                                                                   //(21) 定格暖房能力運転時のヒートポンプサイクルの理論効率 [-]
    /*****************************************************************************************************************************/ 
    alpha_dash_c_hex_C = 0.050 * log((max(V_hs_supply, 400.0) / 3600) / A_f_hex) + 0.073;                                               //(36b) 冷房時の室内熱交換器表面の潜熱伝達率 [kg/(m2・s)] 
    alpha_c_hex_C = alpha_dash_c_hex_C * (c_p_air + c_p_w * 0.010376);                                                                  //(36a) 冷房時の室内熱交換器表面の顕熱伝達率 [W/(m2・K)]
    Theta_sur_f_hex_C = ((Theta_hs_in + Theta_hs_out) / 2) - (c_p_air * rho_air * V_hs_supply * (Theta_hs_in - Theta_hs_out) / 3600) 
                       / (A_e_hex * alpha_c_hex_C);                                                                                     //(32)ndl;    
    Theta_ref_evp_C = max(Theta_sur_f_hex_C, -50.0);                                                                                    //(28) 冷房時の冷媒の蒸発温度 [℃]
    Theta_ref_cnd_C = max(min(max(Theta + 27.4 - 1.35 * Theta_ref_evp_C, Theta), 65.0), Theta_ref_evp_C + 5.0);                         //(27) 冷房時の冷媒の凝縮温度 [℃]
    Theta_ref_SC_C = max(0.772 * Theta_ref_cnd_C - 25.6, 0.0);                                                                          //(29) 冷房時の冷媒の過冷却度 [℃]
    Theta_ref_SH_C = max(0.194 * Theta_ref_cnd_C - 3.86, 0.0);                                                                          //(30) 冷房時の冷媒の過熱度 [℃]
    e_dash_th_C = calc_e_ref_H_th(Theta_ref_evp_C, Theta_ref_cnd_C, Theta_ref_SC_C, Theta_ref_SH_C);                                    //# 4_8_a (1) ヒートポンプサイクルの理論暖房効率 [-];
    e_th_C = e_dash_th_C - 1;                                                                                                           //(21) 定格暖房能力運転時のヒートポンプサイクルの理論効率 [-]
    /*****************************************************************************************************************************/ 
    e_hs_rtd_C = q_hs_rtd_C / (P_hs_rtd_C - P_fan_rtd_C);                                                                               //(12) 定格冷房能力運転時の熱源機の効率 [-]
    e_r_rtd_C = e_hs_rtd_C / e_th_rtd_C;
    e_r_rtd_C = min(max(e_r_rtd_C, 0.0), 1.0);
    e_r_min_C = e_r_rtd_C * 0.65;                                                                                                       //(16) 最小冷房能力運転時のヒートポンプサイクルの理論効率に対する熱源機の効率の比 [-]
    e_r_mid_C = min(max(e_r_rtd_C * 0.95, 0.0), 1.0);                                                                                                     //(14) 中間冷房能力運転時のヒートポンプサイクルの理論効率に対する熱源機の効率の比 [-]
    //e_hs_mid_C = q_hs_mid_C / (P_hs_mid_C - P_fan_mid_C)   #定格能力試験と中間能力試験の値を入力する場合はこちら
    //e_r_mid_C = e_hs_mid_C / e_th_mid_C
    //(10) 冷房時のヒートポンプサイクルの理論効率に対する熱源機の効率の比 [-]

    if(q_hs_C <= q_hs_min_C)                            e_r_C = e_r_min_C - (q_hs_min_C - q_hs_C) * (e_r_min_C / q_hs_min_C);
    else if((q_hs_min_C < q_hs_C) && (q_hs_C <= q_hs_mid_C)) e_r_C = e_r_mid_C - (q_hs_mid_C - q_hs_C) * ((e_r_mid_C - e_r_min_C) / (q_hs_mid_C - q_hs_min_C));
    else if((q_hs_mid_C < q_hs_C) && (q_hs_C <= q_hs_rtd_C)) e_r_C = e_r_rtd_C - (q_hs_rtd_C - q_hs_C) * ((e_r_rtd_C - e_r_mid_C) / (q_hs_rtd_C - q_hs_mid_C));
    else if((q_hs_rtd_C < q_hs_C) && (e_r_rtd_C > 0.4))      e_r_C = max(e_r_rtd_C - (q_hs_C - q_hs_rtd_C) * (e_r_rtd_C / q_hs_rtd_C), 0.4);
    else if((q_hs_rtd_C < q_hs_C) && (e_r_rtd_C <= 0.4))     e_r_C = e_r_rtd_C;
    else                                                     e_r_C = 0.0;

    e_hs_C = e_th_C * e_r_C;                                                                                                            //(8) 冷房時の熱源機の効率 [-]
    E_E_comp_C = ((q_hs_C > 0) && (e_hs_C > 0.0)) ? (q_hs_C / e_hs_C) * 1e-3 : 0.0;                                                                         //(6) 冷房時の圧縮機の消費電力量 [kWh/h]

    double Theta_out = Theta_hs_in - q_hs_CS / (c_p_air * rho_air * V_hs_supply * 3600 / 1e+6);
    double X_out     = X_hs_in     - q_hs_CL / (L_wtr   * rho_air * V_hs_supply * 1e+3 * 3600 / 1e+6);

    return make_tuple(Theta_out, X_out, E_E_comp_C, E_E_fan_C);
}
