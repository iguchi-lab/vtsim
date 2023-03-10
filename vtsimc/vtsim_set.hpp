/***************************************************************************
*  計算定数
***************************************************************************/
#define SOLVE_LU    0
#define SOLVE_SOR   1

#define STEP_P      1e-6        //偏微分時の圧力変化
#define VENT_ERR    1e-6        //換気回路網の許容残差
#define STEP_T      1e-6        //偏微分時の温度変化
#define THRM_ERR    1e-6        //熱回路網の許容残差
#define CONV_ERR    1e-6        //収束許容誤差
#define SOR_RATIO   0.9         //SOR法の緩和係数
#define SOR_ERR     1e-6        //SOR法の許容残差

/***************************************************************************
*  物理定数
***************************************************************************/
#define AIR_CP      1006        //空気の熱容量
#define G           9.81        //重力加速度

//空気密度の計算
double get_rho(double sita){return 353.25 / (sita + 273.15);}
#define RHO20   get_rho(20.0)   //20℃の空気の密度

double  c_p_air        = 1006.0,                //空気の比熱       [J/Kg・K]
        rho_air        = 1.2,                   //空気の密度       [kg/m3]
        L_wtr          = 2500.8 - 2.3668 * 27,  //水の蒸発潜熱     [kJ/kg]
        c_p_w          = 1.846,                 //水蒸気の定圧比熱 [J/Kg・K]
        F              = 101325,                //大気圧          [Pa]
        A_f_hex        = 0.23559,               //室内機熱交換器の前面面積のうち熱交換に有効な面積 [m2]
        A_e_hex        = 6.396;                 //室内機熱交換器の表面積のうち熱交換に有効な面積 [m2] 

/***************************************************************************
*  ノードのフラグ
***************************************************************************/
#define SN_NONE   0           //計算しない
#define SN_CALC   1           //計算する
#define SN_FIX    2           //固定値（計算には利用するが、更新しない）
#define SN_DLY    3           //遅延（熱容量計算用）

/***************************************************************************
*  換気回路網のフラグ
***************************************************************************/
#define VN_SIMPLE   0           //換気回路網：単純開口
#define VN_GAP      1           //換気回路網：隙間
#define VN_FIX      2           //換気回路網：風量固定
#define VN_AIRCON   3           //換気回路網：エアコン=風量固定、換気による熱移動=0
#define VN_FAN      4           //換気回路網：送風ファン、PQ特性

/***************************************************************************
*  熱回路網のフラグ
***************************************************************************/
#define TN_SIMPLE   0           //熱回路網：単純熱回路
#define TN_AIRCON   1           //熱回路網：エアコン、熱量収支付け替え
#define TN_SOLAR    2           //熱回路網：日射取得
#define TN_GROUND   3           //熱回路網：地盤
#define TN_HEATER   4           //熱回路網：発熱

/***************************************************************************
*  エアコンのフラグ
***************************************************************************/
#define AC_OFF      0
#define AC_ON       1

#define AC_AUTO     0           //エアコン：自動
#define AC_HEATING  1           //エアコン：暖房
#define AC_COOLING  2           //エアコン：冷房
#define AC_STOP     3           //エアコン：停止

#define AC_NONE     0
#define AC_RAC      1
#define AC_DUCT_C   2
#define AC_RAC2     3
#define AC_CRIEPI   4