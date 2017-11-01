#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "flint/fmpz_vec.h"
#include "flint/fmpz_poly.h"

//LD_LIBRARY_PATH=/usr/local/lib
//export LD_LIBRARY_PATH

double PI = 3.1415926;
int bit_num;
int sigma;//离散高斯分布标准差
int B;//[-B, B] 近似离散高斯分布
int d;//同态加密f = x^d + 1
int lamda;//系统安全参数
int Max = 10000;//以概率ro 输出x的方法：随机生成 【0， 10000】 之间的数z，如果z小于（小于等于）ro*10000则输出x，否则不输出
int delta;//系统安全参数
int l;
fmpz_t q;
fmpz_t neg_Delta;
fmpz_t t;
fmpz_t Delta_2;
fmpz_t Delta;
fmpz_t T;

//可达到128bit安全
void System_Param()
{
  sigma = 8;
  B = 10*sigma;
  d = 8192;
  bit_num = 128;
  char *str1 = "100000000000000000000000000000000";
  char *str2 = "8000";
  fmpz_set_str(q, str1, 16);//q = 34028236692093846346337460743176821145665536 = pow(2, 128)
  fmpz_set_str(t, str2, 16);//t = pow(2, 15)=32768
  fmpz_cdiv_q(Delta, q, t);
  fmpz_neg(neg_Delta, Delta);
  fmpz_cdiv_q_ui(Delta_2, Delta, 2);  
  fmpz_sqrt(T, q);//pow(2, 64) = 18446744073709551616
  l = fmpz_flog(q, T);
}

//D_{Z, \sigma}
long SampleZ(int seed)
{
  long xx;
  long z;
  double ro;
  long zz;
  int b;

  srand(time(NULL)+seed);

  do
  { 
    xx = rand()%B;
    z = rand()%Max;
    ro = exp(-(PI*xx*xx)/(sigma*sigma));
  }while(!(z < ro*Max));

    b = rand()%2;
    if(b == 0) b = -1;
    xx = xx*b;

    return xx; 
}

//生成（近似）离散高斯分布
void SampleD(fmpz_poly_t v, int seed)
{
  int i; 
  long x;

  for(i = d-1; i >= 0; i--)
     {
        x = SampleZ(i+seed);
       	fmpz_poly_set_coeff_si(v, i, x);
     }     

}

//s和u都是R2上的多项式
void Gen_R2(fmpz_poly_t v, int seed)
{
   int i;      
   srand(time(NULL)+seed);
   int b;
   for(i = d-1; i >= 0; i--)
   {
      b = rand()%2;
      fmpz_poly_set_coeff_ui(v, i, b);
   }

}

void Gen_Rq_div_2(fmpz_t t, int seed)
{
   int i; 
   char *str = (char *)malloc(sizeof(char)*(bit_num-1));
   int b;
   srand(time(NULL)+seed);
   
   for(i = bit_num - 2; i >= 0; i--) str[i] = rand()%2+48;
   fmpz_set_str(t, str, 2);
   b = rand()%2;

   if(b == 0)  b = -1;
   fmpz_mul_si(t, t, b); 

}
/*
void Gen_Rq(fmpz_poly_t v)
{
  int i; 
  fmpz_t x;
  fmpz_init(x);

  flint_rand_t state;
  flint_randinit(state);

  for(i = d-1; i >= 0; i--)
     { 
        fmpz_randtest_mod_signed(x, state, q);
        fmpz_poly_set_coeff_fmpz(v, i, x);
     }    

  fmpz_clear(x);
  flint_randclear(state);

}
*/

//这个的问题是好像没有小一点的数。。计算上应该效率要差一些，但是这个正常的吗？不是很清楚
void Gen_Rq(fmpz_poly_t v)
{
    int i;
    fmpz_t t;
    fmpz_init(t);
    for(i = d-1; i >= 0; i--)
       {
          Gen_Rq_div_2(t, i);
          fmpz_poly_set_coeff_fmpz(v, i, t);

       }

    fmpz_clear(t);
}

void secret_key_Gen(fmpz_poly_t sk)
{ 
   Gen_R2(sk, 1);
}

void public_key_Gen(fmpz_poly_t sk, fmpz_poly_t pk_p0, fmpz_poly_t pk_p1)
{  
    //pk = ( [-(as+e)]_q, a ), a \in R_q, s \in R_2,  e \in D_{Z, sigma}^d
    fmpz_poly_t e;    
    fmpz_poly_t temp;
    fmpz_poly_init(e);
    fmpz_poly_init(temp);

    SampleD(e, 1);
    Gen_Rq(pk_p1);
      
    fmpz_poly_mul(temp, pk_p1, sk);   
    fmpz_poly_add(temp, temp, e); 
    fmpz_poly_neg(temp, temp);
    fmpz_poly_scalar_smod_fmpz(temp, temp, q);
    fmpz_poly_set(pk_p0, temp);
    fmpz_poly_clear(e);
    fmpz_poly_clear(temp);

}

void SH_Keygen(fmpz_poly_t sk, fmpz_poly_t pk_p0, fmpz_poly_t pk_p1, fmpz_poly_t rlk_r0, fmpz_poly_t rlk_r1)
{
   secret_key_Gen(sk);//生成私钥
   public_key_Gen(sk, pk_p0, pk_p1);//生成公钥
   ReKey_Gen(sk, rlk_r0, rlk_r1);//生成Relinear key
        
}

void Encode(int m, fmpz_poly_t pm)
{ 
   int sgn;
   if(m < 0) sgn = -1;
   else sgn = 1;
   int num;
   m = abs(m);
   double bitnum = log2(m);
   int i; 
   int bit; 
   if(bitnum == (int)bitnum) num  = (int)bitnum + 1;
      else num = (int)bitnum + 1; 

   for(i = 0; i < num; i++) 
     {
        bit = m % 2;
        fmpz_poly_set_coeff_ui(pm, i, bit);
        m = m/2;
     }
   fmpz_poly_scalar_mul_si(pm, pm, sgn);
    
}

void SH_Encrypt(fmpz_poly_t m, fmpz_poly_t pk_p0, fmpz_poly_t pk_p1, fmpz_poly_t c0, fmpz_poly_t c1)
{
   fmpz_poly_t u;
   fmpz_poly_t e1;
   fmpz_poly_t e2;
   fmpz_poly_t temp;
   fmpz_poly_t temp1;
   fmpz_poly_init(u);
   fmpz_poly_init(e1);
   fmpz_poly_init(e2);
   fmpz_poly_init(temp);
   fmpz_poly_init(temp1);

   Gen_R2(u, 2);
 
   SampleD(e1, 2);
   SampleD(e2, 3);

   fmpz_poly_mul(temp, pk_p0, u); 

   fmpz_poly_add(temp, temp, e1);
   
   fmpz_poly_scalar_mul_fmpz(temp1, m, Delta);  
  
   fmpz_poly_add(temp, temp, temp1);
 
   fmpz_poly_scalar_smod_fmpz(temp, temp, q);

   fmpz_poly_set(c0, temp);
   
   fmpz_poly_mul(temp, pk_p1, u);

   fmpz_poly_add(temp, temp, e2);
 
   fmpz_poly_scalar_smod_fmpz(c1, temp, q);
   
   fmpz_poly_clear(u);
   fmpz_poly_clear(e1);
   fmpz_poly_clear(e2);
   fmpz_poly_clear(temp);
   fmpz_poly_clear(temp1);

}

void fmpz_poly_nearest_fmpz(fmpz_poly_t v)
{
   long deg = fmpz_poly_degree(v);
   int i;   
   fmpz_t temp;
   fmpz_t temp1;
   fmpz_init(temp);
   fmpz_init(temp1);

   for(i = 0; i <= deg; i++)
      {
         fmpz_poly_get_coeff_fmpz(temp, v, i);

         fmpz_abs(temp1, temp);

         fmpz_fdiv_r(temp1, temp1, Delta);

         if(fmpz_equal(temp, Delta) == 1) {fmpz_set_ui(temp, 1);}
         else if(fmpz_equal(temp, neg_Delta) == 1) {fmpz_set_si(temp, -1);}
         else if(fmpz_is_zero(temp) == 1) {fmpz_set_ui(temp, 0);}

         else if(fmpz_sgn(temp) == -1 && fmpz_cmp(temp1, Delta_2) > 0)  {fmpz_fdiv_q(temp, temp, Delta);}
         else if(fmpz_sgn(temp) == -1 && fmpz_cmp(temp1, Delta_2) <= 0)  {fmpz_cdiv_q(temp, temp, Delta);}
         else if(fmpz_sgn(temp) == 1 && fmpz_cmp(temp1, Delta_2) >= 0)  {fmpz_cdiv_q(temp, temp, Delta);}
         else{fmpz_fdiv_q(temp, temp, Delta);}

         fmpz_poly_set_coeff_fmpz(v, i, temp);           
      }

   fmpz_clear(temp);
   fmpz_clear(temp1); 
}

//c0 = (c_00, c_01), c1 = (c_10, c_11), c = c0+c1

void SH_Add(fmpz_poly_t c_00, fmpz_poly_t c_01, fmpz_poly_t c_10, fmpz_poly_t c_11, fmpz_poly_t c0, fmpz_poly_t c1)
{
    fmpz_poly_t temp;
    fmpz_poly_init(temp);
    fmpz_poly_add(temp, c_00, c_10);
    fmpz_poly_scalar_smod_fmpz(c0, temp, q);   
    
    fmpz_poly_add(temp, c_01, c_11);
    fmpz_poly_scalar_smod_fmpz(c1, temp, q);
    fmpz_poly_clear(temp);

}

void SH_Mul(fmpz_poly_t c_00, fmpz_poly_t c_01, fmpz_poly_t c_10, fmpz_poly_t c_11, fmpz_poly_t c0, fmpz_poly_t c1, fmpz_poly_t c2)
{
   fmpz_poly_t temp1;
   fmpz_poly_t temp2;

   fmpz_poly_init(temp1);
   fmpz_poly_init(temp2);

   fmpz_poly_mul(c0, c_00, c_10);   
   fmpz_poly_nearest_fmpz(c0);

   fmpz_poly_mul(temp1, c_00, c_11);
   fmpz_poly_mul(temp2, c_01, c_10);
   fmpz_poly_add(c1, temp1, temp2);
   fmpz_poly_nearest_fmpz(c1);

   fmpz_poly_mul(c2, c_01, c_11);
   fmpz_poly_nearest_fmpz(c2);

   fmpz_poly_clear(temp1);
   fmpz_poly_clear(temp2);
   
}

void SH_DecMul(fmpz_poly_t sk, fmpz_poly_t c0, fmpz_poly_t c1, fmpz_poly_t c2)
{
  fmpz_poly_t temp;
  fmpz_poly_t temp2;
  fmpz_poly_init(temp);
  fmpz_poly_init(temp2);

  fmpz_poly_mul(temp, c1, sk);

  fmpz_poly_mul(temp2, c2, sk);
  fmpz_poly_mul(temp2, temp2, sk);

  fmpz_poly_add(temp,c0, temp);
  fmpz_poly_add(temp, temp, temp2);

  fmpz_poly_scalar_smod_fmpz(temp, temp, q);   

  fmpz_poly_nearest_fmpz(temp);

  fmpz_poly_scalar_smod_fmpz(temp, temp, t);
   
  fmpz_poly_clear(temp);
  fmpz_poly_clear(temp2);

}

void ReKey_Gen(fmpz_poly_t sk, fmpz_poly_t* rlk_r0, fmpz_poly_t* rlk_r1)
{
   fmpz_poly_t e;
   int i;
   fmpz_poly_t temp1;
   fmpz_t temp2;
   fmpz_poly_t temp3;
   fmpz_poly_init(e);
   fmpz_poly_init(temp1);
   fmpz_init(temp2);
   fmpz_poly_init(temp3);
   
  
   for(i = 0; i < l; i++)
      {
         Gen_Rq(rlk_r1[i]);
         SampleD(e, i+9);
         fmpz_poly_mul(temp1, rlk_r1[i], sk);
         fmpz_poly_add(temp1, temp1, e);
         fmpz_poly_neg(temp1, temp1);
         fmpz_pow_ui(temp2, T, i);
         fmpz_poly_mul(temp3, sk, sk);     
         fmpz_poly_scalar_mul_fmpz(temp3, temp3, temp2);
         fmpz_poly_add(temp1, temp1, temp3);
         fmpz_poly_scalar_smod_fmpz(rlk_r0[i], temp1, q);  
         
      }
   
   fmpz_poly_clear(temp1);
   fmpz_clear(temp2);
   fmpz_poly_clear(temp3);
   fmpz_poly_clear(e);

}
 
//T=sqrt(q), 所以到了i = 2和以后更大的i就没有必要再算了，因为mod q都是0
void  Relinear(fmpz_poly_t* rlk_r0, fmpz_poly_t* rlk_r1, fmpz_poly_t c0, fmpz_poly_t c1, fmpz_poly_t c2, fmpz_poly_t _c0, fmpz_poly_t _c1)
{
   int i; 
   int deg = fmpz_poly_degree(c2);

   fmpz_poly_t *c2_i;
   c2_i = (fmpz_poly_t *)malloc(sizeof(fmpz_poly_t)*l);

   fmpz_t temp;  
   fmpz_t temp1;
   fmpz_poly_t temp2;
   fmpz_poly_t temp3;
 
   for(i = 0; i < l; i++)
      fmpz_poly_init(c2_i[i]);

   fmpz_init(temp);
   fmpz_init(temp1);
   fmpz_poly_init(temp2);
   fmpz_poly_init(temp3);

   for(i = 0; i <= deg; i++)
      {
         fmpz_poly_get_coeff_fmpz(temp, c2, i);
         
         fmpz_abs(temp1, temp);
         if(fmpz_cmp(temp1, T) < 0)
            fmpz_poly_set_coeff_fmpz(c2_i[0], i, temp); 
         else 
            {
               fmpz_fdiv_qr(temp, temp1, temp, T);
               fmpz_poly_set_coeff_fmpz(c2_i[0], i, temp1);
               fmpz_poly_set_coeff_fmpz(c2_i[1], i, temp);
            }

      }

   fmpz_poly_mul(temp2, rlk_r0[0], c2_i[0]);
   fmpz_poly_mul(temp3, rlk_r0[1], c2_i[1]);
  
   fmpz_poly_add(temp2, temp2, temp3);
   fmpz_poly_add(temp2, temp2, c0);
   fmpz_poly_scalar_smod_fmpz(_c0, temp2, q); 

   fmpz_poly_mul(temp2, rlk_r1[0], c2_i[0]);
   fmpz_poly_mul(temp3, rlk_r1[1], c2_i[1]);
  
   fmpz_poly_add(temp2, temp2, temp3);
   fmpz_poly_add(temp2, temp2, c1);
   fmpz_poly_scalar_smod_fmpz(_c1, temp2, q); 

   for(i = 0; i < l; i++)
      fmpz_poly_clear(c2_i[i]);
 
   fmpz_clear(temp);
   fmpz_clear(temp1);
   fmpz_poly_clear(temp2);
   fmpz_poly_clear(temp3);

}

void SH_Decrypt(fmpz_poly_t sk, fmpz_poly_t c0, fmpz_poly_t c1)
{
   fmpz_poly_t temp;
   fmpz_poly_init(temp);  
   fmpz_poly_mul(temp, c1, sk);
   fmpz_poly_add(temp, temp, c0);
   fmpz_poly_scalar_smod_fmpz(temp, temp, q);   
   fmpz_poly_nearest_fmpz(temp);
   fmpz_poly_scalar_smod_fmpz(temp, temp, t);

   printf("\ntemp*t/q is \n");
   fmpz_poly_print(temp);
 
   fmpz_poly_clear(temp);

}

struct timespec diff(struct timespec start,struct timespec end)
{
    struct timespec temp;
    if((end.tv_nsec-start.tv_nsec)<0)
    {
        temp.tv_sec = end.tv_sec-start.tv_sec-1;
        temp.tv_nsec = 1000000000 + end.tv_nsec-start.tv_nsec;
    }
    else
    {
        temp.tv_sec = end.tv_sec-start.tv_sec;
        temp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }
    return temp;
}

int main()
{
  fmpz_poly_t pk_p0;
  fmpz_poly_t pk_p1;
  fmpz_poly_t sk;
  fmpz_poly_t* rlk_r0;
  fmpz_poly_t* rlk_r1;
  int m;
  int i;
  fmpz_poly_t m0;
  fmpz_poly_t m1;
  fmpz_poly_t c_00;
  fmpz_poly_t c_01;
  fmpz_poly_t c_10;
  fmpz_poly_t c_11;
  fmpz_poly_t c0;
  fmpz_poly_t c1;
  fmpz_poly_t c2;
  fmpz_poly_t _c0;
  fmpz_poly_t _c1;
  struct timespec tpstart;
  struct timespec tpend;
  struct timespec tp;
  long timedif; 
  
  fmpz_init(q);
  fmpz_init(neg_Delta);
  fmpz_init(Delta_2);
  fmpz_init(t);
  fmpz_init(Delta);
  fmpz_init(T);

  fmpz_poly_init(pk_p0);
  fmpz_poly_init(pk_p1);
  fmpz_poly_init(sk);
  fmpz_poly_init(m0);
  fmpz_poly_init(m1);
  fmpz_poly_init(c_00);
  fmpz_poly_init(c_01);
  fmpz_poly_init(c_10);
  fmpz_poly_init(c_11);
  fmpz_poly_init(c0);
  fmpz_poly_init(c1);
  fmpz_poly_init(c2); 
  fmpz_poly_init(_c0);
  fmpz_poly_init(_c1);

  System_Param();

  rlk_r0 = (fmpz_poly_t *)malloc(sizeof(fmpz_poly_t)*(l));  
  rlk_r1 = (fmpz_poly_t *)malloc(sizeof(fmpz_poly_t)*(l)); 

  for(i = 0; i < l; i++)
     {
        fmpz_poly_init(rlk_r0[i]);  
        fmpz_poly_init(rlk_r1[i]);
     }

  clock_gettime(CLOCK_MONOTONIC, &tpstart);
  SH_Keygen(sk, pk_p0, pk_p1, rlk_r0, rlk_r1);
  clock_gettime(CLOCK_MONOTONIC, &tpend);
  tp = diff(tpstart,tpend);
  printf("\nKEYGEN%ds:%ldns \n", tp.tv_sec, tp.tv_nsec);

  scanf("%d", &m);
  Encode(m, m0);
  scanf("%d", &m);
  Encode(m, m1);

  clock_gettime(CLOCK_MONOTONIC, &tpstart);
  SH_Encrypt(m0, pk_p0, pk_p1, c_00, c_01);
  clock_gettime(CLOCK_MONOTONIC, &tpend);
  tp = diff(tpstart,tpend);
  printf("\nENC%ds:%ldns \n", tp.tv_sec, tp.tv_nsec);

  SH_Encrypt(m1, pk_p0, pk_p1, c_10, c_11);

  clock_gettime(CLOCK_MONOTONIC, &tpstart);
  SH_Decrypt(sk, c_00, c_01);
  clock_gettime(CLOCK_MONOTONIC, &tpend);
  tp = diff(tpstart,tpend);
  printf("\nDEC%ds:%ldns \n", tp.tv_sec, tp.tv_nsec);

  SH_Decrypt(sk, c_10, c_11); 
 
  clock_gettime(CLOCK_MONOTONIC, &tpstart);
  SH_Add(c_00, c_01, c_10, c_11, c0, c1);
  clock_gettime(CLOCK_MONOTONIC, &tpend);
  tp = diff(tpstart,tpend);
  printf("\nADD%ds:%ldns \n", tp.tv_sec, tp.tv_nsec);

  SH_Decrypt(sk, c0, c1);
   
  clock_gettime(CLOCK_MONOTONIC, &tpstart);
  SH_Mul(c_00, c_01, c_10, c_11, c0, c1, c2);
  clock_gettime(CLOCK_MONOTONIC, &tpend);
  tp = diff(tpstart,tpend);
  printf("\nMUL%ds:%ldns \n", tp.tv_sec, tp.tv_nsec);


 // SH_DecMul(sk, c0, c1, c2);
  
  Relinear(rlk_r0, rlk_r1, c0, c1, c2, _c0, _c1);
  SH_Decrypt(sk, _c0, _c1);

   
  fmpz_poly_clear(sk);  
  fmpz_poly_clear(pk_p0);
  fmpz_poly_clear(pk_p1);
  fmpz_clear(q);  
  fmpz_clear(neg_Delta);
  fmpz_clear(Delta_2);  
  fmpz_clear(t);
  fmpz_clear(Delta);
  fmpz_clear(T);
  fmpz_poly_clear(c_00);
  fmpz_poly_clear(c_01);
  fmpz_poly_clear(c_10);
  fmpz_poly_clear(c_11);
  fmpz_poly_clear(c0);
  fmpz_poly_clear(c1);
  fmpz_poly_clear(c2);
  fmpz_poly_clear(_c0);
  fmpz_poly_clear(_c1);
  fmpz_poly_clear(m0);
  fmpz_poly_clear(m1);
  
  for(i = 0; i < l; i++) 
     {
       fmpz_poly_clear(rlk_r0[i]);
       fmpz_poly_clear(rlk_r1[i]);
     }

  return 0;
}
