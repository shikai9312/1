CONTANT
=====
## 模拟模型和条件
* 打开Energy模型和k-e模型RNG
* 空气入口速度5m/s，冷却壁面恒温20℃
* 出口pressure-outlet，回流温度为293.15K
* 壁面convection，h=25W/m^2，无穷远处温度为293.15K
* Solution Method选择couple，二阶迎风，Green-Gauss Cell Based
* 残差设置1e-8

## pathline的设置
scale=0.05，color by Velocity，Release from各个面

## 两相流模型
## UDF编辑
* UDF导入使用complier的模式。必须保存成.c的文件才能导入
* 不知道为什么，这个文件的real变量必须在main程序之外
* Umax需要调整
```
/************************************************************************
UDF that computes the particle breakage PDF
*************************************************************************/
#include "udf.h"
#include "sg_pb.h"
#include "sg_mphase.h"

#define kb 1.381e-23
#define T 343.15
#define u 1.7894e-05
#define v 1.48e-05
#define c 20.99948   /*湍流聚并常数*/ 
#define Umax 4e-18   /*颗粒间排斥势能，需要调整*/
#define den 2100
#define e 40      /*单位质量能量耗散率*/
real agg_kernel;
real a1;   /*布朗聚并*/
real a2;   /*湍流聚并*/
real R0;
real f;   /*捕集效率*/
real mp;   /*聚并后颗粒质量*/

DEFINE_PB_COALESCENCE_RATE(aggregation_kernel,cell,thread,d_1,thread_2,d_2)
{   
	a1 = 2*kb*T*pow((d_1+d_2),2.0)/(3*u*(d_1*d_2));	
	R0 = (d_1+d_2)/4;	
	mp = 0.523599*(pow(d_1,3.0)+pow(d_2,3.0))*den;
	f = (9*Umax*v/(4*R0*R0*e*mp)+1)*exp(-9*Umax*v/(4*R0*R0*e*mp));
	a2 = c*f*pow(R0,3.0)*sqrt(e/v);
	agg_kernel = sqrt( pow(a1,2.0) + pow(a2,2.0));
    return agg_kernel;
}

```
* 或者更精确地，获取当地位置的单位质量能量耗散率来进行计算。（这里直接读取的颗粒1的能量耗散率。因为考虑到在当地1和2的耗散率是相同的。能量耗散率越大，聚并效果越好）

```
/************************************************************************
UDF that computes the particle breakage PDF
*************************************************************************/
#include "udf.h"
#include "sg_pb.h"
#include "sg_mphase.h"

#define kb 1.381e-23
#define T 343.15
#define u 1.7894e-05
#define v 1.48e-05
#define c 20.99948   /*湍流聚并常数*/ 
#define Umax 4e-18   /*颗粒间排斥势能，需要调整*/
#define den 2100
real e;      /*单位质量能量耗散率*/
real agg_kernel;
real a1;   /*布朗聚并*/
real a2;   /*湍流聚并*/
real R0;
real f;   /*捕集效率*/
real mp;   /*聚并后颗粒质量*/

DEFINE_PB_COALESCENCE_RATE(aggregation_kernel,cell,thread,d_1,thread_2,d_2)
{   
	Thread *tm = THREAD_SUPER_THREAD(thread);
	e = C_O(cell,tm);
	a1 = 2*kb*T*pow((d_1+d_2),2.0)/(3*u*(d_1*d_2));	
	R0 = (d_1+d_2)/4;	
	mp = 0.523599*(pow(d_1,3.0)+pow(d_2,3.0))*den;
	f = (9*Umax*v/(4*R0*R0*e*mp)+1)*exp(-9*Umax*v/(4*R0*R0*e*mp));
	a2 = c*f*pow(R0,3.0)*sqrt(e/v);
	agg_kernel = sqrt( pow(a1,2.0) + pow(a2,2.0));
        return agg_kernel;
}

```

