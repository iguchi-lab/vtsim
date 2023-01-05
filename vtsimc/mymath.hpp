#include <math.h>

using namespace std;

tuple<vector<vector<double>>, vector<double>> pivoting(vector<vector<double>> a, vector<double> b, long size){     // ピボット選択
    for(int pivot = 0; pivot < size; pivot++){                                                          // 各列で 一番値が大きい行を 探す
        int     max_row =   pivot;
        double  max_val =   0;
        for (int row = pivot; row < size; row++){
            if (fabs(a[row][pivot]) > max_val){
                max_val =   fabs(a[row][pivot]);                                                        // 一番値が大きい行
                max_row =   row;
            }
        }
        if (max_row != pivot){                                                                          // 一番値が大きい行と入れ替え
            double tmp;
            for (int col = 0; col < size; col++){
                tmp             =   a[max_row][col];
                a[max_row][col] =   a[pivot][col];
                a[pivot][col]   =   tmp;
            }
            tmp         =   b[max_row];
            b[max_row]  =   b[pivot];
            b[pivot]    =   tmp;
        }
    }
    return make_tuple(a, b);
}

//参考HP    http://fornext1119.web.fc2.com/NumericOperation/vol_06/Text/10_05_14.xhtml
vector<double> LU(vector<vector<double>> a, vector<double> b, long size){
    vector<double> x(size, 0.0), y(size, 0.0);
    tie(a, b) = pivoting(a, b, size);
    //LOG_PRINT("a :");
    for (int pivot = 0; pivot < size - 1; pivot++){
        for (int row = pivot + 1; row < size; row++){
            double s = a[row][pivot] / a[pivot][pivot];
            
            for (int col = pivot; col < size; col++){ 
                a[row][col] -= a[pivot][col] * s;          // これが 上三角行列
                //LOG_CONTENTS("c_a[" << row << "][" << col << "] = " << a[row][col] <<  ", ");
            }
            a[row][pivot] = s;                                                                  // これが 下三角行列
            //LOG_CONTENTS("p_a[" << row << "][" << pivot << "] = " << a[row][pivot] <<  ", ");
        }
    }
    //LOG_CONTENTS(endl);

    for (int row = 0; row < size; row++){
        for (int col = 0; col < row; col++)         b[row] -= a[row][col] * y[col];
        y[row] = b[row];
    }

    for (int row = size - 1; row >= 0; row--){
        for (int col = size - 1; col > row; col--)  y[row] -= a[row][col] * x[col];
        x[row] = y[row] / a[row][row];
    }
    return x;
}
     
vector<double> SOR(vector<vector<double>> a, vector<double> b, long size, double sor_ratio, double sor_err){     //SOR法による計算
    double dd, d2, sum_a;
    vector<double> x(size, 0.0);                                    //xの初期化（出力用）
    tie(a, b) = pivoting(a, b, size);
    do{
        d2 = 0;                                                     //残差の初期化
        for(int i = 0; i < size; i++){
            sum_a = b[i];                                           //sumaをbで初期化
            for(int j = 0; j < size; j++)  
                if(i != j)      sum_a -= a[i][j] * x[j];
            dd = sor_ratio * (sum_a / a[i][i] - x[i]);              //dpの増分の計算
            x[i] += dd;                                             //dpの更新
            d2 += sqrt(pow(dd, 2.0));                               //増分から残差を計算
        }
    }while(sor_err < d2);                             //残差が小さくなればループを抜ける
    return x;
}