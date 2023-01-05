using namespace std;

class Node{
public:
    int i, s_i;
    tuple<int, int, int, int> flag;                         //換気、濃度、熱計算フラグ                
    vector<double> p, c, t, x;                              //圧力、濃度、温度、絶対湿度
    vector<double> m, mx;                                   //粉塵発生量、発湿量
    vector<double> h_sr, h_inp;                             //日射量、発熱量
    vector<double> v;                                       //気積     
    vector<double> beta;                                    //沈着率

    Node(long length, int i_, tuple<int, int, int, int> flag_){
        i = i_;
        flag = flag_;
        p.assign(length, 0.0);
        c.assign(length, 0.0);
        t.assign(length, 20.0);
        x.assign(length, 0.0);
    }

};