
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>


using namespace std;

vector<double> KIR_linear(const double& a,const double& t,const double& h,const vector<double>& w_old,const int& iteration){
    vector<double> w_new (w_old.size(),0);
    if(a>0){
    for(unsigned long int i = iteration;i<w_new.size();++i) {
        w_new[i] = w_old[i] + (abs(a) - a) / 2 * t / h * (w_old[i + 1] - w_old[i]) -
                   (a + abs(a)) / 2 * t / h * (w_old[i] - w_old[i - 1]);
    }}
    else{
            for (unsigned long int i = 0; i < w_new.size()-iteration;++i){
                w_new[i] = w_old[i] + (abs(a) - a) / 2 * t / h * (w_old[i + 1] - w_old[i]) -
                           (a + abs(a)) / 2 * t / h * (w_old[i] - w_old[i - 1]);
            }
        }
    return w_new;
}

vector<double> KIR_nonlinear(const double& a,const double& t,const double& h,const vector<double>& w_old,const int& iteration){
    vector<double> w_new (w_old.size(),0);
    if(a>0){
        for(unsigned long int i = iteration;i<w_new.size();++i) {
            w_new[i] = w_old[i] + (abs(a) - a) / 2 * t / h * (w_old[i + 1] - w_old[i]) -
                       (a + abs(a)) / 2 * t / h * (w_old[i] - w_old[i - 1])+abs(a)/a*t;;
        }}
    else{
        for (unsigned long int i = 0; i < w_new.size()-iteration;++i){
            w_new[i] = w_old[i] + (abs(a) - a) / 2 * t / h * (w_old[i + 1] - w_old[i]) -
                       (a + abs(a)) / 2 * t / h * (w_old[i] - w_old[i - 1])+abs(a)/a*t;;
        }
    }
    return w_new;
}

int main() {
    //вводим время в котором хотим посмотреть решение
    //без ГУ ис такой матрицей оно может быть не больше чем пи/2
    //я удлинил для удобства промежуток
    double t_interest;
    //ЧТобы не хранить большое количество данных программа дает решение в конкретный момент времени
    cout<<"Задайте момент времени"<<endl;
    cin>>t_interest;
    if(t_interest>M_PI/2) {
        cout<<"вышли за треугольник"<<endl;
        exit(0);
    }
    //отрезок
    double L = M_PI;
    //шаг по пространству
    double h = 0.01;
    //шаг по времени
    double t = 0.5*h;
    //создаем сетки
    double nh = L/h;
    vector<double> hi (nh,0);

    for (int i =0;i<nh;++i){
        hi[i] = i*h;
    }
    //вектора W
    vector<double> w1(nh,0);
    vector<double> w2(nh,0);
    vector<double> w3(nh,0);
    vector<double> w4(nh,0);
    //начальные условия на сразу два случая азадчи
    for(int i =0;i<nh;++i){
        w1[i] = (sin(i*h)+cos(i*h));
        w2[i] = sin(i*h)-cos(i*h);
        w3[i] = (sin(i*h)+cos(i*h));
        w4[i] = -sin(i*h)+cos(i*h);
    }
    //применяем разные схемы КИР
    for(int i =0;i<t_interest/t;++i){
        w1 = KIR_linear(1,t,h,w1,i+1);
        w2 = KIR_linear(-1,t,h,w2,i+1);
        w3 = KIR_nonlinear(1,t,h,w3,i+1);
        w4 = KIR_nonlinear(-1,t,h,w4,i+1);
    }
    //Программа берет от вас номер точки массива данных в указанный момент времени и выдает значение функций в точке n*h
    //где h это шаг

    //однородная система
    while(true){
        double n;
        cout<<"Однородная система.Задайте номер точки (координата n*h).Для выхода -1"<<endl;
        cin>>n;
        if((n*h<0)||(n*h>L)){
            cout<<"выход за границы уравнения"<<endl;
            break;
        }
        if((w1[n]==0)||(w2[n]==0)){
            cout<<"выход за треуголльник"<<endl;

        }
        double u_calc = (w1[n]+w2[n])/2;
        double v_calc = (w1[n]-w2[n])/2;
        double u_anal = sin(n*h)*(cos(t_interest)+sin(t_interest));
        double v_anal = cos(h*n)*(cos(t_interest)-sin(t_interest));
        cout<<"u_calc = "<<u_calc<<"||||| u_anal = "<<u_anal<<endl;
        cout<<"v_calc = "<<v_calc<<"||||| v_anal = "<<v_anal<<endl;
    }

    //неоднородная система
    while(true){
        double n;
        cout<<"Неоднородная система.Задайте номер точки (координата n*h).Для выхода -1"<<endl;
        cin>>n;
        if((n*h<0)||(n*h>L)){
            cout<<"выход за границы уравнения"<<endl;
            break;
        }
        if((w3[n]==0)||(w4[n]==0)){
            cout<<"выход за треуголльник"<<endl;

        }
        double u_calc = (w3[n]-w4[n])/2;
        double v_calc = (w3[n]+w4[n])/2;
        double u_anal = sin(n*h)*(cos(t_interest)+sin(t_interest))+t_interest;
        double v_anal = cos(h*n)*(cos(t_interest)-sin(t_interest));
        cout<<"u_calc = "<<u_calc<<"||||| u_anal = "<<u_anal<<endl;
        cout<<"v_calc = "<<v_calc<<"||||| v_anal = "<<v_anal<<endl;
    }

    //запишем данные для графика

    string path_zero ="data_zero.txt" ;
    string path_one ="data_one.txt" ;
    ofstream out_zero,out_one;          // поток для записи
    out_zero.open(path_zero);// окрываем файл для записи
    if (out_zero.is_open())
    {
        cout<<"open data_zero correctly"<<endl;
        for(int i =0;i<w1.size();++i){
            out_zero<<(w1[i]+w2[i])/2<<" "<<sin(i*h)*(cos(t_interest)+sin(t_interest))<<" "<<(w1[i]-w2[i])/2<<" "<<cos(h*i)*(cos(t_interest)-sin(t_interest))<<" "<<i*h<<endl;
        }
        out_zero<<t_interest;
        out_zero.close();

    }
    out_one.open(path_one);
    if (out_one.is_open())
    {
        cout<<"open data_one correctly"<<endl;
        for(int i =0;i<w1.size();++i){
            out_one<<(w3[i]-w4[i])/2<<" "<<sin(i*h)*(cos(t_interest)+sin(t_interest))+t_interest<<" "<<(w3[i]+w4[i])/2<<" "<<cos(h*i)*(cos(t_interest)-sin(t_interest))<<" "<<i*h<<endl;
        }
        out_one<<t_interest;
        out_one.close();


    }



    return 0;
}
