#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

float p0= 0.5;
float x0 = 1.0;
double x_min = -2.0;
double x_max = 2.0;
double p_min = -2.0;
double p_max = 2.0;
const int n_x = 100;
const int n_p = 100;
double epsilon = 1.0;
double d_t = 0.01;


double x_values[n_x];
double p_values[n_p];
double force[n_x];
double main_func[n_x][n_p];
double func_one[n_x][n_p];
double func_two[n_x][n_p];


//do the function of % but with negative numbers
int b_c(int value)
{
    if (value < 0)
    {
        value = n_x + value;
    }
    return value;
}


// ---------------------------- waterbag array-----------------------------------------------------
void water_bag(double function[n_x][n_p])
{
    for (int j = 0; j < n_p; j++)
    {
        for (int i = 0; i < n_x; i++)
        {
            if (abs(x_values[i]) < x0)
            {
                if (abs(p_values[j]) < p0)
                {
                    function[i][j] = 1.0 / (4.0 * x0 * p0);
                } else {
                    function[i][j] = 0.0;
                }
            }
        }
    } 
}
// ----------------------------Creation of x and p values array-----------------------------------------------------
void init_xvalues()
{
    for (int i = 0; i < n_x; i++)
    {
        x_values[i] = x_min + (x_max - x_min) * i / (n_x - 1);
    }
}
void init_pvalues()
{
    for (int i = 0; i < n_p; i++)
    {
        p_values[i] = p_min + ((p_max - p_min) * i) / (n_p - 1);
    }
}

//--------------------------------FUNCIONES DEL CUBIC SPLINE-------------------------------------------------
double Alpha(double x, int i)
{
    return (x_values[(i + 1) % n_x] - x) / (x_values[(i + 1) % n_x] - x_values[i]);
}
double Beta(double x, int i)
{
    return 1 - Alpha(x, i);
}
double Gamma(double x, int i)
{
double alpha = Alpha(x, i);
    return ((pow(alpha, 3) - alpha) / 6) * pow((x_values[(i + 1) % n_x] - x_values[i]), 2);
}
double Delta(double x, int i)
{
    double beta = Beta(x, i);
    return ((pow(beta, 3) - beta) / 6) * pow((x_values[(i + 1) % n_x] - x_values[i]), 2);
}

//---------------------------------SECONDS DERIVATES----------------------------------------------------------

//second derivative in x
double sec_drv_x(double func_array[n_x][n_p],int i,int j)
{
    return -(1/560)*func_array[b_c(i-4)][j] + (8/315)*func_array[b_c(i-3)][j] - (1/5)*func_array[b_c(i-2)][j] + (8/5)*func_array[b_c(i-1)][j] - (205/72)*func_array[(i)%n_x][j] + (8/5)*func_array[(i+1)%n_x][j] - (1/5)*func_array[(i+2)%n_x][j] + (8/315)*func_array[(i+3)%n_x][j] - (1/560)*func_array[(i+4)%n_x][j];
}
//second derivative in x with i+1
double sec_drv_x_plus(double func_array[n_x][n_p],int i,int j)
{
    return -(1/560)*func_array[b_c(i+1-4)][j] + (8/315)*func_array[b_c(i+1-3)][j] - (1/5)*func_array[b_c(i+1-2)][j] + (8/5)*func_array[i+1-1][j] - (205/72)*func_array[(i+1)%n_x][j] + (8/5)*func_array[(i+1+1)%n_x][j] - (1/5)*func_array[(i+1+2)%n_x][j] + (8/315)*func_array[(i+1+3)%n_x][j] - (1/560)*func_array[(i+1+4)%n_x][j];
}

//second derivative in p
double sec_drv_p(double func_array[n_x][n_p],int i,int j)
{
    return -(1/560)*func_array[i][b_c(j-4)] + (8/315)*func_array[i][b_c(j-3)] - (1/5)*func_array[i][b_c(j-2)] + (8/5)*func_array[i][b_c(j-1)] - (205/72)*func_array[i][(j)%n_x] + (8/5)*func_array[i][(j+1)%n_x] - (1/5)*func_array[i][(j+2)%n_x] + (8/315)*func_array[i][(j+3)%n_x] - (1/560)*func_array[i][(j+4)%n_x];
}
//second derivative in p with j+1
double sec_drv_p_plus(double func_array[n_x][n_p],int i,int j)
{
    return -(1/560)*func_array[i][b_c(j+1-4)] + (8/315)*func_array[i][b_c(j+1-3)] - (1/5)*func_array[i][b_c(j+1-2)] + (8/5)*func_array[i][j+1-1] - (205/72)*func_array[i][(j+1)%n_x] + (8/5)*func_array[i][(j+1+1)%n_x] - (1/5)*func_array[i][(j+1+2)%n_x] + (8/315)*func_array[i][(j+1+3)%n_x] - (1/560)*func_array[i][(j+1+4)%n_x];
}



//--------------------------------FUNCIONES DEL CUBIC SPLINE-------------------------------------------------

double func_csp_x(double func_array[n_x][n_p],double x,int i,int j)
{
    return Alpha(x,i)*func_array[i][j]+ Beta(x,i)*func_array[(i+1)%n_x][j]+ Gamma(x,i)*sec_drv_x(func_array,i,j) + Delta(x,i)*sec_drv_x_plus(func_array,i,j)   ;
}
double func_csp_p(double func_array[n_x][n_p],double p,int i,int j)
{
    return  Alpha(p,j)*func_array[i][j] + Beta(p,j)*func_array[i][(j+1)%n_p] + Gamma(p,j)*sec_drv_p(func_array,i,j) + Delta(p,j)*sec_drv_p_plus(func_array,i,j);

}

//------------DEBE SER UNA INTEGRACION NUMERICA, CASO RAPIDO, TRAPECIOS, CASO LENTO, CUADRATURA--------------

//trapesoid integral k=1 to n-1, in this case 1 to n_x-1.
double sum_fxj(double fab[n_x][n_p],int x) 
{
    double sum = 0;
    for (int j = 1; j < n_p-1; j++)
    {
        sum = sum + fab[x][j];
    }
    return sum;
}
//Is the rightest therms of the double trapesoid integral to the space integral with "b" variable x_values[k]
double sum_f_sum_b(double fab[n_x][n_p],int k) // k do the function of n-1 in the loop
{
    double sum = 0;
    for (int i = 1; i < k; i++)
    {
        for (int j = 1; j < n_p-1; j++)
        {
            sum = sum + (((fab[i][0]+fab[i][n_p-1])*(1/2)) + sum_fxj(fab,i));
        }
    }
    return sum;
}
//Is the rightest therms of the double trapesoid integral to the space integral with "a" variable x_values[k]
double sum_f_sum_a(double fab[n_x][n_p],int k) // k do the function of n-1 in the loop
{
    double sum = 0;
    for (int i = k; i < n_x-1; i++)
    {
        for (int j = 1; j < n_p-1; j++)
        {
            sum = sum + (((fab[i][0]+fab[i][n_p-1])*(1/2)) + sum_fxj(fab,i));
        }
    }
    return sum;
}

//Double trapezoid integral and add result to force array
void integration(double force[n_x],double f[n_x][n_p])
{
    for (int k=0 ; k < n_x ; k++)
    {
        double hx_b = (x_values[k]-x_min)/ (k+1);
        double hx_a = (x_max-x_values[k])/ (n_x-k);
        double hp = (p_max - p_min)/ n_p ;

        if ( k < 2)
        {
            force[k] = hp*(epsilon)*(hx_b*(-1)*((1/4)*(f[0][0]+f[0][n_p-1]+f[k][0]+f[k][n_p-1])) ) + (hx_a*( (1/4)*(f[k][0]+f[k][n_p-1] + f[n_x-1][k]+f[k][n_p-1]) ));
        }
        else
        {
            force[k] = (hx_b*hp*( ((1/4)*( f[0][0]+f[0][n_p-1] + 2*sum_fxj(f,0) + f[k][0]+f[k][n_p-1] + sum_fxj(f,k) )) +sum_f_sum_b(f,k)  )) + ( hx_a*hp*( (1/4)*((f[k][0]+f[k][n_p-1]+2*sum_fxj(f,k)) + (f[n_x-1][0]+f[n_x-1][n_p-1]+2*sum_fxj(f,n_x-1))) + sum_f_sum_a(f,k))) ;
        }
    }
}

//-----------------------------------FUNCIONES DEL LOOP--------------------------------------------------------

void advection_x(double func_array[n_x][n_p], double func_one[n_x][n_p])
{
    for (int j = 0; j < n_p; j++) {
        for (int i = 0; i < n_x; i++) {
            func_one[i][j] = func_csp_x(func_array,x_values[i]-((p_values[j]*d_t)/2),i,j);
        }
    }
}
void advection_p(double func_array[n_x][n_p], double func_one[n_x][n_p],double force[n_x])
{
    for (int j = 0; j < n_p; j++) {
        for (int i = 0; i < n_x; i++) {
            func_one[i][j] = func_csp_p(func_array,p_values[j] + (force[i] * d_t),i,j);
        }
    }
}

//-----------------------------------

void reset(double func[n_x])
{
    for (int i = 0; i < n_p; i++) {
            func[i] = 0.0;
    }
}



//-----------------------------------MAIN--------------------------------------------------------


int main()
{ 

    init_xvalues();//inicia los valores de x // SE RECOMIENDA CAMBIAR POR ALGUN CONSTRUCTOR Y NO INICIALIZAR
    init_pvalues();//inicia los valores de p  
    water_bag(main_func);// crea la funcion inicial
        
    // printf("main_func \n");
    //  for (int i = 0; i < n_x; i++) {
    //     for (int j = 0; j < n_p; j++) {
    //         //cout << main_func[(i+1)%n_x][j];
    //         printf("%f ", main_func[i][j]);
    //     }
    //     printf("\n");
    // } 

    cout<< main_func[n_x/2][n_p/2]<<endl;
    int loop=0;
    for(int t = 0; t < 40; t++)
    {
        advection_x(main_func,func_one);// realiza el primer paso, first advetion in x.
        reset(force);// resetea el array force
        integration(force,func_one);// realiza la integracion y la guarda en force

        advection_p(func_one,func_two,force);// realiza el segundo paso, second advetion in p.

        advection_x(func_two,main_func);// realiza el tercer paso, third advetion in x.
        cout<<"func_one"<<func_one[n_x/2][n_p/2]<<endl;
        cout<<"func_two"<<func_two[n_x/2][n_p/2]<<endl;
        cout<<"main_func"<<main_func[n_x/2][n_p/2]<<endl;
        loop++;
        cout<<"---------------loop "<<loop<<endl;
    }
  



    // printf("main after loop %d \n" ,loop);
    //  for (int i = 0; i < n_x; i++) {
    //     for (int j = 0; j < n_p; j++) {
    //         //cout << main_func[(i+1)%n_x][j];
    //         printf("%f ", main_func[i][j]);
    //     }
    //     printf("\n");
    // }   
    // cout<<"v_dxmaxima"<<(x_values[2]-x_values[1])/d_t <<endl;
    // for (int i = 0; i < n_x; i++)
    // {

    //     cout<<force[i]<<endl;
    // }
    //cout<<force[4]<<endl;

    
    return 0;
}