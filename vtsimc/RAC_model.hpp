
double  table_3[5][6] =   {{-0.00236,  0.01324,  0.08418, -0.47143, -1.16944,   6.54886},
                           { 0.00427, -0.02392, -0.19226,  0.94213,  2.58632, -12.85618},
                           {-0.00275,  0.01542,  0.14947, -0.68303, -2.03594,  10.60561},
                           { 0.00063, -0.00351, -0.02865,  0.10522,  0.37336,  -1.09499},
                           {-0.00005,  0.00028,  0.00184, -0.01090, -0.09609,   0.59229}},
        table_4_A[5][3] = {{-0.000056,  0.000786,  0.071625},
                           {-0.000145,  0.003337, -0.143643},
                           {-0.000240, -0.029471,  1.954343},
                           {-0.000035, -0.050909,  1.389751},
                           { 0.0,       0.0,       0.076800}},
        table_4_B[5][3] = {{ 0.000108, -0.035658,  3.063873},
                           {-0.000017,  0.062546, -5.471556},
                           {-0.000245, -0.025126,  4.057590},
                           { 0.000323, -0.021166,  0.575459},
                           { 0.0,       0.000330,  0.047500}},
        table_4_C[5][3] = {{-0.001465, -0.030500,  1.920431},
                           { 0.002824,  0.041081, -1.835302},
                           {-0.001929, -0.009738,  1.582898},
                           { 0.000616, -0.014239,  0.546204},
                           { 0.0,      -0.000110,  0.023100}},
        table_5[5][6] =   {{0.00000, 0.00000,  0.00000,  0.00000,  0.00000,  0.00000},
                           {0.00000, 0.00000, -0.00036,  0.05080, -0.20346,  0.47765},
                           {0.00000, 0.00000,  0.00227, -0.03952,  0.04115,  0.23099},
                           {0.00000, 0.00000, -0.00911,  0.07102,  0.14950, -1.07335},
                           {0.00000, 0.00000,  0.00044, -0.00214, -0.06250,  0.35150}},
        table_6_A[5][3] = {{-0.0004078,  0.01035,  -0.03248},
                           { 0.0,        0.04099,  -0.818809},
                           { 0.0,       -0.04615,   2.10666},
                           { 0.0013382, -0.01179,  -0.41778},
                           { 0.0000000, -0.00102,   0.09270}},
        table_6_B[5][3] = {{-0.000056,  -0.003539, -0.430566},
                           { 0.0,        0.015237,  1.188850},
                           { 0.0,        0.000527, -0.304645},
                           {-0.000179,   0.020543,  0.130373},
                           { 0.0,        0.000240,  0.013500}},
        table_6_C[5][3] = {{-0.0001598,  0.004848,  0.047097},
                           { 0.0,        0.016675,  0.362141},
                           { 0.0,       -0.008134, -0.023535},
                           {-0.0000772,  0.012558,  0.056185},
                           { 0.0,       -0.000110,  0.010300}};

double  C_af_H = 0.8,       //室内機吹き出し風量に関する暖房出力補正係数 [-]
        C_df_H = 1.0,       //デフロストに関する暖房出力補正係数 [-]
        C_af_C = 0.85,      //室内機吹き出し風量に関する冷房出力補正係数 [-] 
        C_hm_C = 1.15,      //室内機吸い込み湿度に関する冷房出力補正係数 [-]
        SHF_L_min_c = 0.4;  //処理冷房全熱負荷に対する処理冷房顕熱負荷の最小値 [-]

/*********************************************************************************************
 暖房
**********************************************************************************************/
tuple<double, double, double> calc_a_eq3(double q_r_max_H, double q_rtd_C){
    double a0, a1, a2, b0, b1, b2, c0, c1, c2;
    b2 =  0.000181 * q_rtd_C * 1e-3 - 0.000184;
    b1 =  0.002322 * q_rtd_C * 1e-3 + 0.013904;
    b0 =  0.003556 * q_rtd_C * 1e-3 + 0.993431;
    c2 = -0.000173 * q_rtd_C * 1e-3 + 0.000367;
    c1 = -0.003980 * q_rtd_C * 1e-3 + 0.003983;
    c0 = -0.002870 * q_rtd_C * 1e-3 + 0.006376;
    a2 = b2 * q_r_max_H + c2;
    a1 = b1 * q_r_max_H + c1;
    a0 = b0 * q_r_max_H + c0;
    return make_tuple(a2, a1, a0);  //(3a),(3b),(3c)
}

double calc_Q_r_max_H(double q_rtd_C, double q_r_max_H, double Theta_ex_d_t){
    double a0, a1, a2;
    tie(a2, a1, a0) = calc_a_eq3(q_r_max_H, q_rtd_C);
    return a2 * pow(Theta_ex_d_t - 7, 2) + a1 * (Theta_ex_d_t - 7) + a0;    //(2)
}

double calc_p_i_eq8(int i, double q_rtd_C){
    double s_i, t_i;
    s_i = table_3[int(4 - floor(i / 10))][(2 - (i % 10)) * 2];
    t_i = table_3[int(4 - floor(i / 10))][(2 - (i % 10)) * 2 + 1];
    return s_i * q_rtd_C * 1e-3 + t_i;  //(8)
}

double calc_p_i_eq9(int i, double q_rtd_C){
    double p_i_A, p_i_B, p_i_C;
    p_i_A = table_4_A[int(4 - floor(i / 10))][(2 - (i % 10))];
    p_i_B = table_4_B[int(4 - floor(i / 10))][(2 - (i % 10))];
    p_i_C = table_4_C[int(4 - floor(i / 10))][(2 - (i % 10))];
    if(q_rtd_C <= 2200)                       
        return p_i_A;
    else if((2200 < q_rtd_C) && (q_rtd_C <= 4000))    
        return p_i_A * ((4000 - q_rtd_C) / (4000 - 2200)) + p_i_B * ((q_rtd_C - 2200) / (4000 - 2200));
    else if((4000 < q_rtd_C) && (q_rtd_C < 7100))     
        return p_i_B * ((7100 - q_rtd_C) / (7100 - 4000)) + p_i_C * ((q_rtd_C - 4000) / (7100 - 4000));
    else if(7100 <= q_rtd_C)                       
        return p_i_C;
    return 0.0;     //(9a),(9b),(9c),(9d)
}

tuple<double, double, double, double, double>calc_a_eq7(double Theta_ex, double q_rtd_C, bool dualcompressor){
    double a0, a1, a2, a3, a4, p00, p01, p02, p10, p11, p12, p20, p21, p22, p30, p31, p32, p40, p41, p42;  
    if(dualcompressor == false){
        p42 = calc_p_i_eq8(42, q_rtd_C);
        p41 = calc_p_i_eq8(41, q_rtd_C);
        p40 = calc_p_i_eq8(40, q_rtd_C);
        p32 = calc_p_i_eq8(32, q_rtd_C);
        p31 = calc_p_i_eq8(31, q_rtd_C);
        p30 = calc_p_i_eq8(30, q_rtd_C);
        p22 = calc_p_i_eq8(22, q_rtd_C);
        p21 = calc_p_i_eq8(21, q_rtd_C);
        p20 = calc_p_i_eq8(20, q_rtd_C);
        p12 = calc_p_i_eq8(12, q_rtd_C);
        p11 = calc_p_i_eq8(11, q_rtd_C);
        p10 = calc_p_i_eq8(10, q_rtd_C);
        p02 = calc_p_i_eq8(2,  q_rtd_C);
        p01 = calc_p_i_eq8(1,  q_rtd_C);
        p00 = calc_p_i_eq8(0,  q_rtd_C);
    }
    else if(dualcompressor == true){
        p42 = calc_p_i_eq9(42, q_rtd_C);
        p41 = calc_p_i_eq9(41, q_rtd_C);
        p40 = calc_p_i_eq9(40, q_rtd_C);
        p32 = calc_p_i_eq9(32, q_rtd_C);
        p31 = calc_p_i_eq9(31, q_rtd_C);
        p30 = calc_p_i_eq9(30, q_rtd_C);
        p22 = calc_p_i_eq9(22, q_rtd_C);
        p21 = calc_p_i_eq9(21, q_rtd_C);
        p20 = calc_p_i_eq9(20, q_rtd_C);
        p12 = calc_p_i_eq9(12, q_rtd_C);
        p11 = calc_p_i_eq9(11, q_rtd_C);
        p10 = calc_p_i_eq9(10, q_rtd_C);
        p02 = calc_p_i_eq9(2,  q_rtd_C);
        p01 = calc_p_i_eq9(1,  q_rtd_C);
        p00 = calc_p_i_eq9(0,  q_rtd_C);
    }
    a4 = p42 * pow(Theta_ex, 2) + p41 * Theta_ex + p40 * 1;
    a3 = p32 * pow(Theta_ex, 2) + p31 * Theta_ex + p30 * 1;
    a2 = p22 * pow(Theta_ex, 2) + p21 * Theta_ex + p20 * 1;
    a1 = p12 * pow(Theta_ex, 2) + p11 * Theta_ex + p10 * 1;
    a0 = p02 * pow(Theta_ex, 2) + p01 * Theta_ex + p00 * 1;
    return make_tuple(a0, a1, a2, a3, a4);  //(7)
}

double calc_f_H_Theta(double x, double Theta_ex, double q_rtd_C, bool dualcompressor){
    double a0, a1, a2, a3, a4;
    tie(a0, a1, a2, a3, a4) = calc_a_eq7(Theta_ex, q_rtd_C, dualcompressor);
    return a4 * pow(x, 4) + a3 * pow(x, 3) + a2 * pow(x, 2) + a1 * x + a0; //(6)
}

/*********************************************************************************************
 冷房
**********************************************************************************************/
tuple<double, double, double> calc_a_eq13(double q_r_max_C, double q_rtd_C){
    double a0, a1, a2, b0, b1, b2, c0, c1, c2;
    b2 =  0.000812 * q_rtd_C * 1e-3 - 0.001480;
    b1 =  0.003527 * q_rtd_C * 1e-3 - 0.023000;
    b0 = -0.011490 * q_rtd_C * 1e-3 + 1.024328;
    c2 = -0.000350 * q_rtd_C * 1e-3 + 0.000800;
    c1 = -0.001280 * q_rtd_C * 1e-3 + 0.003621;
    c0 =  0.004772 * q_rtd_C * 1e-3 - 0.011170;
    a2 = b2 * q_r_max_C + c2;
    a1 = b1 * q_r_max_C + c1;
    a0 = b0 * q_r_max_C + c0;
    return make_tuple(a2, a1, a0); //(13)
}

double calc_Q_r_max_C(double q_r_max_C, double q_rtd_C, double Theta_ex_d_t){
    double a2, a1, a0;
    tie(a2, a1, a0) = calc_a_eq13(q_r_max_C, q_rtd_C);
    return a2 * pow(Theta_ex_d_t - 35, 2) + a1 * (Theta_ex_d_t - 35) + a0;  //(12)
}

double calc_p_i_eq23(int i, double q_rtd_C){
    double s_i, t_i;
    s_i = table_5[int(4 - floor(i / 10))][(2 - (i % 10)) * 2];
    t_i = table_5[int(4 - floor(i / 10))][(2 - (i % 10)) * 2 + 1];
    return s_i * q_rtd_C * 1e-3 + t_i;  //(23)
}

double calc_p_i_eq24(int i, double q_rtd_C){
    double p_i_A, p_i_B, p_i_C;
    p_i_A = table_6_A[int(4 - floor(i / 10))][(2 - (i % 10))];
    p_i_B = table_6_B[int(4 - floor(i / 10))][(2 - (i % 10))];
    p_i_C = table_6_C[int(4 - floor(i / 10))][(2 - (i % 10))];
    if(q_rtd_C <= 2200)                    
        return p_i_A;
    else if((2200 < q_rtd_C) && (q_rtd_C <= 4000)) 
        return p_i_A * ((4000 - q_rtd_C) / (4000 - 2200)) + p_i_B * ((q_rtd_C - 2200) / (4000 - 2200));
    else if((4000 <= q_rtd_C) && (q_rtd_C < 7100)) 
        return p_i_B * ((7100 - q_rtd_C) / (7100 - 4000)) + p_i_C * ((q_rtd_C - 4000) / (7100 - 4000));
    else if(7100 <= q_rtd_C)                    
        return p_i_C;
    return 0.0; //(24)
}

tuple<double, double, double, double , double> calc_a_eq22(double Theta_ex, double q_rtd_C, bool dualcompressor){
    double a0, a1, a2, a3, a4, p00, p01, p02, p10, p11, p12, p20, p21, p22, p30, p31, p32, p40, p41, p42;  
    if(dualcompressor == false){
        p42 = calc_p_i_eq23(42, q_rtd_C);
        p41 = calc_p_i_eq23(41, q_rtd_C);
        p40 = calc_p_i_eq23(40, q_rtd_C);
        p32 = calc_p_i_eq23(32, q_rtd_C);
        p31 = calc_p_i_eq23(31, q_rtd_C);
        p30 = calc_p_i_eq23(30, q_rtd_C);
        p22 = calc_p_i_eq23(22, q_rtd_C);
        p21 = calc_p_i_eq23(21, q_rtd_C);
        p20 = calc_p_i_eq23(20, q_rtd_C);
        p12 = calc_p_i_eq23(12, q_rtd_C);
        p11 = calc_p_i_eq23(11, q_rtd_C);
        p10 = calc_p_i_eq23(10, q_rtd_C);
        p02 = calc_p_i_eq23(2,  q_rtd_C);
        p01 = calc_p_i_eq23(1,  q_rtd_C);
        p00 = calc_p_i_eq23(0,  q_rtd_C);   
    }
    else if (dualcompressor == true){
        p42 = calc_p_i_eq24(42, q_rtd_C);
        p41 = calc_p_i_eq24(41, q_rtd_C);
        p40 = calc_p_i_eq24(40, q_rtd_C);
        p32 = calc_p_i_eq24(32, q_rtd_C);
        p31 = calc_p_i_eq24(31, q_rtd_C);
        p30 = calc_p_i_eq24(30, q_rtd_C);
        p22 = calc_p_i_eq24(22, q_rtd_C);
        p21 = calc_p_i_eq24(21, q_rtd_C);
        p20 = calc_p_i_eq24(20, q_rtd_C);
        p12 = calc_p_i_eq24(12, q_rtd_C);
        p11 = calc_p_i_eq24(11, q_rtd_C);
        p10 = calc_p_i_eq24(10, q_rtd_C);
        p02 = calc_p_i_eq24(2,  q_rtd_C);
        p01 = calc_p_i_eq24(1,  q_rtd_C);
        p00 = calc_p_i_eq24(0,  q_rtd_C);   
    }
    a4 = p42 * pow(Theta_ex, 2) + p41 * Theta_ex + p40 * 1;
    a3 = p32 * pow(Theta_ex, 2) + p31 * Theta_ex + p30 * 1;
    a2 = p22 * pow(Theta_ex, 2) + p21 * Theta_ex + p20 * 1;
    a1 = p12 * pow(Theta_ex, 2) + p11 * Theta_ex + p10 * 1;
    a0 = p02 * pow(Theta_ex, 2) + p01 * Theta_ex + p00 * 1;
    return make_tuple(a0, a1, a2, a3, a4);  //(22)
}

double calc_f_C_Theta(double x, double Theta_ex, double q_rtd_C, bool dualcompressor){
    double a0, a1, a2, a3, a4;
    tie(a0, a1, a2, a3, a4) = calc_a_eq22(Theta_ex, q_rtd_C, dualcompressor);
    return a4 * pow(x,  4) + a3 * pow(x, 3) + a2 * pow(x, 2) + a1 * x + a0; //(21)
}