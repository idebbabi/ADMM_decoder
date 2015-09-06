


void sort15_rank_order_reg_modif(double llr[], double rllr[], int ipos[], int rpos[])
{
    const float x0  = llr[0];  const float x1  = llr[1];
    const float x2  = llr[2];  const float x3  = llr[3];
    const float x4  = llr[4];  const float x5  = llr[5];
    const float x6  = llr[6];  const float x7  = llr[7];
    const float x8  = llr[8];  const float x9  = llr[9];
    const float x10 = llr[10]; const float x11 = llr[11];
    const float x12 = llr[12]; const float x13 = llr[13];
    const float x14 = llr[14]; 

    const int o0 = (x0> x1)+(x0 >x2)+(x0 >x3)+(x0 >x4)+(x0 >x5)+(x0 >x6)+(x0 >x7)+(x0 >x8)+(x0 >x9)+(x0 >x10)+(x0 >x11)+(x0 >x12)+(x0 >x13)+ (x0 >x14) ;
    const int o1 = (x1>=x0)+(x1 >x2)+(x1 >x3)+(x1 >x4)+(x1 >x5)+(x1 >x6)+(x1 >x7)+(x1 >x8)+(x1 >x9)+(x1 >x10)+(x1 >x11)+(x1 >x12)+(x1 >x13)+ (x1 >x14) ;
    const int o2 = (x2>=x0)+(x2>=x1)+(x2 >x3)+(x2 >x4)+(x2 >x5)+(x2 >x6)+(x2 >x7)+(x2 >x8)+(x2 >x9)+(x2 >x10)+(x2 >x11)+(x2 >x12)+(x2 >x13)+ (x2 >x14) ;
    const int o3 = (x3>=x0)+(x3>=x1)+(x3>=x2)+(x3 >x4)+(x3 >x5)+(x3 >x6)+(x3 >x7)+(x3 >x8)+(x3 >x9)+(x3 >x10)+(x3 >x11)+(x3 >x12)+(x3 >x13)+ (x3 >x14) ;
    const int o4 = (x4>=x0)+(x4>=x1)+(x4>=x2)+(x4>=x3)+(x4 >x5)+(x4 >x6)+(x4 >x7)+(x4 >x8)+(x4 >x9)+(x4 >x10)+(x4 >x11)+(x4 >x12)+(x4 >x13)+ (x4 >x14) ;
    const int o5 = (x5>=x0)+(x5>=x1)+(x5>=x2)+(x5>=x3)+(x5>=x4)+(x5 >x6)+(x5 >x7)+(x5 >x8)+(x5 >x9)+(x5 >x10)+(x5 >x11)+(x5 >x12)+(x5 >x13)+ (x5 >x14) ;
    const int o6 = (x6>=x0)+(x6>=x1)+(x6>=x2)+(x6>=x3)+(x6>=x4)+(x6>=x5)+(x6 >x7)+(x6 >x8)+(x6 >x9)+(x6 >x10)+(x6 >x11)+(x6 >x12)+(x6 >x13)+ (x6 >x14) ;

    const int o7 = (x7 >=x0)+(x7 >=x1)+(x7 >=x2)+(x7 >=x3)+(x7 >=x4)+(x7 >=x5)+(x7  >=x6)+(x7  > x8)+(x7 >  x9)+(x7  > x10)+(x7  > x11)+(x7 > x12)+(x7  > x13) + (x7  >x14) ;
    const int o8 = (x8 >=x0)+(x8 >=x1)+(x8 >=x2)+(x8 >=x3)+(x8 >=x4)+(x8 >=x5)+(x8  >=x6)+(x8  >=x7)+(x8 >  x9)+(x8  > x10)+(x8  > x11)+(x8 > x12)+(x8  > x13) + (x8  >x14) ;
    const int o9 = (x9 >=x0)+(x9 >=x1)+(x9 >=x2)+(x9 >=x3)+(x9 >=x4)+(x9 >=x5)+(x9  >=x6)+(x9  >=x7)+(x9 >= x8)+(x9  > x10)+(x9  > x11)+(x9 > x12)+(x9  > x13) + (x9  >x14) ;
    const int o10= (x10>=x0)+(x10>=x1)+(x10>=x2)+(x10>=x3)+(x10>=x4)+(x10>=x5)+(x10 >=x6)+(x10 >=x7)+(x10>= x8)+(x10 >= x9)+(x10 > x11)+(x10> x12)+(x10 > x13) + (x10 >x14) ;
    const int o11= (x11>=x0)+(x11>=x1)+(x11>=x2)+(x11>=x3)+(x11>=x4)+(x11>=x5)+(x11 >=x6)+(x11 >=x7)+(x11>= x8)+(x11 >= x9)+(x11 >=x10)+(x11> x12)+(x11 > x13) + (x11 >x14) ;
    const int o12= (x12>=x0)+(x12>=x1)+(x12>=x2)+(x12>=x3)+(x12>=x4)+(x12>=x5)+(x12 >=x6)+(x12 >=x7)+(x12>= x8)+(x12 >= x9)+(x12 >=x10)+(x12>=x11)+(x12 > x13) + (x12 >x14) ;
    const int o13= (x13>=x0)+(x13>=x1)+(x13>=x2)+(x13>=x3)+(x13>=x4)+(x13>=x5)+(x13 >=x6)+(x13 >=x7)+(x13>= x8)+(x13 >= x9)+(x13 >=x10)+(x13>=x11)+(x13 >=x12) + (x12 >x14) ;

    const int o14 = 105 - (o0 + o1 + o2 + o3 + o4 + o5 + o6 +o7 +o8 +o9 +o10 +o11 +o12+ o13);

    rllr[o0]=x0; rllr[o1]=x1; rllr[o2]=x2; rllr[o3]=x3;
    rllr[o4]=x4; rllr[o5]=x5; rllr[o6]=x6; rllr[o7]=x7;
    rllr[o8]=x8; rllr[o9]=x9; rllr[o10]=x10; rllr[o11]=x11;
    rllr[o12]=x12; rllr[o13]=x13; rllr[o14]=x14;


    rpos[o0]=ipos[0]; rpos[o1]=ipos[1]; rpos[o2]=ipos[2]; rpos[o3]=ipos[3];
    rpos[o4]=ipos[4]; rpos[o5]=ipos[5]; rpos[o6]=ipos[6]; rpos[o7]=ipos[7];
    rpos[o8]=ipos[8]; rpos[o9]=ipos[9]; rpos[o10]=ipos[10]; rpos[o11]=ipos[11];
    rpos[o12]=ipos[12]; rpos[o13]=ipos[13]; rpos[o14]=ipos[14];

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void sort15_rank_order_reg_modif(float llr[], float rllr[], int ipos[], int rpos[])
{
    const float x0  = llr[0];  const float x1  = llr[1];
    const float x2  = llr[2];  const float x3  = llr[3];
    const float x4  = llr[4];  const float x5  = llr[5];
    const float x6  = llr[6];  const float x7  = llr[7];
    const float x8  = llr[8];  const float x9  = llr[9];
    const float x10 = llr[10]; const float x11 = llr[11];
    const float x12 = llr[12]; const float x13 = llr[13];
    const float x14 = llr[14]; 

    const int o0 = (x0> x1)+(x0 >x2)+(x0 >x3)+(x0 >x4)+(x0 >x5)+(x0 >x6)+(x0 >x7)+(x0 >x8)+(x0 >x9)+(x0 >x10)+(x0 >x11)+(x0 >x12)+(x0 >x13)+ (x0 >x14) ;
    const int o1 = (x1>=x0)+(x1 >x2)+(x1 >x3)+(x1 >x4)+(x1 >x5)+(x1 >x6)+(x1 >x7)+(x1 >x8)+(x1 >x9)+(x1 >x10)+(x1 >x11)+(x1 >x12)+(x1 >x13)+ (x1 >x14) ;
    const int o2 = (x2>=x0)+(x2>=x1)+(x2 >x3)+(x2 >x4)+(x2 >x5)+(x2 >x6)+(x2 >x7)+(x2 >x8)+(x2 >x9)+(x2 >x10)+(x2 >x11)+(x2 >x12)+(x2 >x13)+ (x2 >x14) ;
    const int o3 = (x3>=x0)+(x3>=x1)+(x3>=x2)+(x3 >x4)+(x3 >x5)+(x3 >x6)+(x3 >x7)+(x3 >x8)+(x3 >x9)+(x3 >x10)+(x3 >x11)+(x3 >x12)+(x3 >x13)+ (x3 >x14) ;
    const int o4 = (x4>=x0)+(x4>=x1)+(x4>=x2)+(x4>=x3)+(x4 >x5)+(x4 >x6)+(x4 >x7)+(x4 >x8)+(x4 >x9)+(x4 >x10)+(x4 >x11)+(x4 >x12)+(x4 >x13)+ (x4 >x14) ;
    const int o5 = (x5>=x0)+(x5>=x1)+(x5>=x2)+(x5>=x3)+(x5>=x4)+(x5 >x6)+(x5 >x7)+(x5 >x8)+(x5 >x9)+(x5 >x10)+(x5 >x11)+(x5 >x12)+(x5 >x13)+ (x5 >x14) ;
    const int o6 = (x6>=x0)+(x6>=x1)+(x6>=x2)+(x6>=x3)+(x6>=x4)+(x6>=x5)+(x6 >x7)+(x6 >x8)+(x6 >x9)+(x6 >x10)+(x6 >x11)+(x6 >x12)+(x6 >x13)+ (x6 >x14) ;

    const int o7 = (x7 >=x0)+(x7 >=x1)+(x7 >=x2)+(x7 >=x3)+(x7 >=x4)+(x7 >=x5)+(x7  >=x6)+(x7  > x8)+(x7 >  x9)+(x7  > x10)+(x7  > x11)+(x7 > x12)+(x7  > x13) + (x7  >x14) ;
    const int o8 = (x8 >=x0)+(x8 >=x1)+(x8 >=x2)+(x8 >=x3)+(x8 >=x4)+(x8 >=x5)+(x8  >=x6)+(x8  >=x7)+(x8 >  x9)+(x8  > x10)+(x8  > x11)+(x8 > x12)+(x8  > x13) + (x8  >x14) ;
    const int o9 = (x9 >=x0)+(x9 >=x1)+(x9 >=x2)+(x9 >=x3)+(x9 >=x4)+(x9 >=x5)+(x9  >=x6)+(x9  >=x7)+(x9 >= x8)+(x9  > x10)+(x9  > x11)+(x9 > x12)+(x9  > x13) + (x9  >x14) ;
    const int o10= (x10>=x0)+(x10>=x1)+(x10>=x2)+(x10>=x3)+(x10>=x4)+(x10>=x5)+(x10 >=x6)+(x10 >=x7)+(x10>= x8)+(x10 >= x9)+(x10 > x11)+(x10> x12)+(x10 > x13) + (x10 >x14) ;
    const int o11= (x11>=x0)+(x11>=x1)+(x11>=x2)+(x11>=x3)+(x11>=x4)+(x11>=x5)+(x11 >=x6)+(x11 >=x7)+(x11>= x8)+(x11 >= x9)+(x11 >=x10)+(x11> x12)+(x11 > x13) + (x11 >x14) ;
    const int o12= (x12>=x0)+(x12>=x1)+(x12>=x2)+(x12>=x3)+(x12>=x4)+(x12>=x5)+(x12 >=x6)+(x12 >=x7)+(x12>= x8)+(x12 >= x9)+(x12 >=x10)+(x12>=x11)+(x12 > x13) + (x12 >x14) ;
    const int o13= (x13>=x0)+(x13>=x1)+(x13>=x2)+(x13>=x3)+(x13>=x4)+(x13>=x5)+(x13 >=x6)+(x13 >=x7)+(x13>= x8)+(x13 >= x9)+(x13 >=x10)+(x13>=x11)+(x13 >=x12) + (x12 >x14) ;

    const int o14 = 105 - (o0 + o1 + o2 + o3 + o4 + o5 + o6 +o7 +o8 +o9 +o10 +o11 +o12+ o13);

    rllr[o0]=x0; rllr[o1]=x1; rllr[o2]=x2; rllr[o3]=x3;
    rllr[o4]=x4; rllr[o5]=x5; rllr[o6]=x6; rllr[o7]=x7;
    rllr[o8]=x8; rllr[o9]=x9; rllr[o10]=x10; rllr[o11]=x11;
    rllr[o12]=x12; rllr[o13]=x13; rllr[o14]=x14;


    rpos[o0]=ipos[0]; rpos[o1]=ipos[1]; rpos[o2]=ipos[2]; rpos[o3]=ipos[3];
    rpos[o4]=ipos[4]; rpos[o5]=ipos[5]; rpos[o6]=ipos[6]; rpos[o7]=ipos[7];
    rpos[o8]=ipos[8]; rpos[o9]=ipos[9]; rpos[o10]=ipos[10]; rpos[o11]=ipos[11];
    rpos[o12]=ipos[12]; rpos[o13]=ipos[13]; rpos[o14]=ipos[14];

}


