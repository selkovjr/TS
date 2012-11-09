/* Copyright (C) 2011 Ion Torrent Systems, Inc. All Rights Reserved */
// This is a collection of pre-defined smoothing parameters for Tikhonov smoothing of the traces
// (APB)
int tik05_frame_cnt = 25;
int tik05_denominator = 16384;
int tik05_i_low[25] = {
       0,      0,      0,      0,      0,      0,      0,      0,      1,      2, 
       3,      4,      5,      6,      7,      8,      9,     10,     11,     12, 
      13,     14,     15,     16,     17
};
int tik05_i_high[25] = {
       7,      8,      9,     10,     11,     12,     13,     14,     15,     16, 
      17,     18,     19,     20,     21,     22,     23,     24,     24,     24, 
      24,     24,     24,     24,     24
};
int tik05_imatrix[25][25] = {
{

   12882,   3967,    154,   -409,   -193,    -34,      9,      8,      0,      0, 
       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

    3967,   9069,   3404,    370,   -249,   -151,    -35,      3,      6,      0, 
       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

     154,   3404,   9285,   3563,    412,   -250,   -156,    -37,      3,      6, 
       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

    -409,    370,   3563,   9327,   3563,    407,   -252,   -157,    -37,      3, 
       6,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

    -193,   -249,    412,   3563,   9321,   3560,    407,   -252,   -157,    -37, 
       3,      6,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

     -34,   -151,   -250,    407,   3560,   9322,   3560,    407,   -252,   -157, 
     -37,      3,      6,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

       9,    -35,   -156,   -252,    407,   3560,   9321,   3560,    407,   -252, 
    -157,    -37,      3,      6,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

       8,      3,    -37,   -157,   -252,    407,   3560,   9322,   3561,    407, 
    -252,   -157,    -37,      3,      6,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

       0,      6,      3,    -37,   -157,   -252,    407,   3560,   9322,   3561, 
     407,   -252,   -157,    -37,      3,      6,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

       0,      0,      6,      3,    -37,   -157,   -252,    407,   3561,   9322, 
    3561,    407,   -252,   -157,    -37,      3,      6,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

       0,      0,      0,      6,      3,    -37,   -157,   -252,    407,   3561, 
    9322,   3561,    407,   -252,   -157,    -37,      3,      6,      0,      0, 
       0,      0,      0,      0,      0
},
{

       0,      0,      0,      0,      6,      3,    -37,   -157,   -252,    407, 
    3561,   9322,   3561,    407,   -252,   -157,    -37,      3,      6,      0, 
       0,      0,      0,      0,      0
},
{

       0,      0,      0,      0,      0,      6,      3,    -37,   -157,   -252, 
     407,   3561,   9322,   3561,    407,   -252,   -157,    -37,      3,      6, 
       0,      0,      0,      0,      0
},
{

       0,      0,      0,      0,      0,      0,      6,      3,    -37,   -157, 
    -252,    407,   3561,   9322,   3561,    407,   -252,   -157,    -37,      3, 
       6,      0,      0,      0,      0
},
{

       0,      0,      0,      0,      0,      0,      0,      6,      3,    -37, 
    -157,   -252,    407,   3561,   9322,   3561,    407,   -252,   -157,    -37, 
       3,      6,      0,      0,      0
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      6,      3, 
     -37,   -157,   -252,    407,   3561,   9322,   3561,    407,   -252,   -157, 
     -37,      3,      6,      0,      0
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      6, 
       3,    -37,   -157,   -252,    407,   3561,   9322,   3560,    407,   -252, 
    -157,    -37,      3,      6,      0
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       6,      3,    -37,   -157,   -252,    407,   3561,   9321,   3560,    407, 
    -252,   -157,    -37,      3,      9
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      6,      3,    -37,   -157,   -252,    407,   3560,   9321,   3561, 
     407,   -253,   -159,    -39,     16
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      6,      3,    -37,   -157,   -252,    407,   3561,   9322, 
    3560,    405,   -258,   -163,    -13
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      6,      3,    -37,   -157,   -252,    407,   3560, 
    9322,   3562,    409,   -255,   -184
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      6,      3,    -37,   -157,   -253,    405, 
    3562,   9339,   3613,    447,   -544
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0,      6,      3,    -37,   -159,   -258, 
     409,   3613,   9489,   3726,   -408
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0,      0,      6,      3,    -39,   -163, 
    -255,    447,   3726,   9574,   3084
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0,      0,      0,      9,     16,    -13, 
    -184,   -544,   -408,   3085,  14424
}
};
int tik075_frame_cnt = 25;
int tik075_denominator = 16384;
int tik075_i_low[25] = {
       0,      0,      0,      0,      0,      0,      0,      0,      1,      2, 
       3,      4,      5,      6,      7,      8,      9,     10,     11,     12, 
      13,     14,     15,     16,     17
};
int tik075_i_high[25] = {
       7,      8,      9,     10,     11,     12,     13,     14,     15,     16, 
      17,     18,     19,     20,     21,     22,     23,     24,     24,     24, 
      24,     24,     24,     24,     24
};
int tik075_imatrix[25][25] = {
{

   11369,   4915,    956,   -304,   -360,   -166,    -36,      8,      0,      0, 
       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

    4916,   7407,   3655,    900,   -110,   -230,   -122,    -33,      2,      0, 
       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

     957,   3655,   7351,   3849,   1030,    -66,   -227,   -127,    -37,      0, 
       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

    -304,    900,   3849,   7480,   3893,   1033,    -71,   -231,   -129,    -37, 
       1,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

    -360,   -110,   1030,   3893,   7484,   3887,   1029,    -73,   -231,   -129, 
     -37,      1,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

    -166,   -230,    -66,   1033,   3887,   7480,   3886,   1029,    -73,   -231, 
    -129,    -37,      1,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

     -36,   -122,   -227,    -71,   1029,   3886,   7479,   3886,   1030,    -73, 
    -231,   -129,    -37,      1,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

       8,    -33,   -127,   -231,    -73,   1029,   3886,   7479,   3888,   1030, 
     -73,   -231,   -129,    -37,      1,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

       0,      2,    -37,   -129,   -231,    -73,   1029,   3886,   7485,   3889, 
    1030,    -73,   -231,   -129,    -37,      1,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

       0,      0,      0,    -37,   -129,   -231,    -73,   1029,   3888,   7485, 
    3889,   1030,    -73,   -231,   -129,    -37,      1,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

       0,      0,      0,      1,    -37,   -129,   -231,    -73,   1030,   3889, 
    7484,   3889,   1030,    -73,   -231,   -129,    -37,      1,      0,      0, 
       0,      0,      0,      0,      0
},
{

       0,      0,      0,      0,      1,    -37,   -129,   -231,    -73,   1030, 
    3889,   7484,   3889,   1030,    -73,   -231,   -129,    -37,      1,      0, 
       0,      0,      0,      0,      0
},
{

       0,      0,      0,      0,      0,      1,    -37,   -129,   -231,    -73, 
    1030,   3889,   7484,   3889,   1030,    -73,   -231,   -129,    -37,      1, 
       0,      0,      0,      0,      0
},
{

       0,      0,      0,      0,      0,      0,      1,    -37,   -129,   -231, 
     -73,   1030,   3889,   7484,   3889,   1030,    -73,   -231,   -129,    -37, 
       1,      0,      0,      0,      0
},
{

       0,      0,      0,      0,      0,      0,      0,      1,    -37,   -129, 
    -231,    -73,   1030,   3889,   7484,   3889,   1030,    -73,   -231,   -129, 
     -37,      1,      0,      0,      0
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      1,    -37, 
    -129,   -231,    -73,   1030,   3889,   7484,   3889,   1029,    -73,   -231, 
    -129,    -37,      1,      0,      0
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      1, 
     -37,   -129,   -231,    -73,   1030,   3889,   7487,   3886,   1029,    -73, 
    -231,   -129,    -38,      1,      0
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       1,    -37,   -129,   -231,    -73,   1030,   3889,   7479,   3886,   1029, 
     -73,   -233,   -133,    -38,     21
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      1,    -37,   -129,   -231,    -73,   1030,   3886,   7479,   3886, 
    1028,    -76,   -238,   -132,     -9
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      1,    -37,   -129,   -231,    -73,   1029,   3886,   7480, 
    3886,   1029,    -74,   -238,   -145
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      1,    -37,   -129,   -232,    -73,   1028,   3886, 
    7486,   3908,   1063,    -79,   -440
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      1,    -37,   -130,   -233,    -76,   1029, 
    3908,   7550,   4010,   1049,   -688
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0,      1,    -38,   -133,   -238,    -74, 
    1063,   4010,   7710,   3988,     95
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0,      0,      1,    -38,   -132,   -238, 
     -79,   1049,   3988,   7714,   4121
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0,      0,      0,     21,     -9,   -145, 
    -439,   -688,     95,   4119,  13429
}
};
int tik10_frame_cnt = 25;
int tik10_denominator = 16384;
int tik10_i_low[25] = {
       0,      0,      0,      0,      0,      0,      0,      0,      1,      2, 
       3,      4,      5,      6,      7,      8,      9,     10,     11,     12, 
      13,     14,     15,     16,     17
};
int tik10_i_high[25] = {
       7,      8,      9,     10,     11,     12,     13,     14,     15,     16, 
      17,     18,     19,     20,     21,     22,     23,     24,     24,     24, 
      24,     24,     24,     24,     24
};
int tik10_imatrix[25][25] = {
{

   10256,   5323,   1628,     -7,   -382,   -285,   -125,    -28,      0,      0, 
       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

    5327,   6552,   3687,   1254,     89,   -222,   -188,    -89,    -24,      0, 
       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

    1629,   3688,   6176,   3784,   1413,    187,   -186,   -183,    -94,    -28, 
       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

      -7,   1254,   3784,   6337,   3881,   1449,    191,   -191,   -188,    -97, 
     -29,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

    -382,     90,   1413,   3881,   6375,   3886,   1444,    186,   -194,   -188, 
     -97,    -29,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

    -285,   -222,    187,   1449,   3886,   6367,   3881,   1442,    186,   -194, 
    -188,    -97,    -29,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

    -126,   -188,   -186,    191,   1444,   3881,   6365,   3880,   1443,    186, 
    -194,   -188,    -97,    -29,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

     -28,    -89,   -183,   -191,    186,   1442,   3880,   6364,   3882,   1443, 
     186,   -194,   -188,    -97,    -29,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

       0,    -24,    -94,   -188,   -194,    186,   1442,   3881,   6370,   3884, 
    1443,    186,   -194,   -188,    -97,    -29,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

       0,      0,    -28,    -97,   -188,   -194,    186,   1442,   3882,   6373, 
    3884,   1443,    186,   -194,   -188,    -97,    -29,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

       0,      0,      0,    -29,    -97,   -188,   -193,    186,   1443,   3884, 
    6374,   3884,   1443,    186,   -194,   -188,    -97,    -29,      0,      0, 
       0,      0,      0,      0,      0
},
{

       0,      0,      0,      0,    -29,    -96,   -188,   -193,    186,   1443, 
    3884,   6374,   3884,   1443,    186,   -194,   -188,    -96,    -29,      0, 
       0,      0,      0,      0,      0
},
{

       0,      0,      0,      0,      0,    -29,    -96,   -188,   -194,    186, 
    1443,   3884,   6374,   3884,   1443,    186,   -194,   -188,    -96,    -29, 
       0,      0,      0,      0,      0
},
{

       0,      0,      0,      0,      0,      0,    -29,    -96,   -188,   -194, 
     186,   1443,   3884,   6374,   3884,   1444,    186,   -193,   -188,    -96, 
     -29,      0,      0,      0,      0
},
{

       0,      0,      0,      0,      0,      0,      0,    -29,    -97,   -188, 
    -194,    186,   1443,   3884,   6374,   3885,   1444,    186,   -193,   -188, 
     -97,    -29,      0,      0,      0
},
{

       0,      0,      0,      0,      0,      0,      0,      0,    -29,    -97, 
    -188,   -194,    186,   1443,   3884,   6373,   3886,   1442,    186,   -194, 
    -189,    -98,    -30,      0,      0
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,    -29, 
     -97,   -188,   -194,    186,   1443,   3885,   6376,   3881,   1442,    186, 
    -195,   -191,   -101,    -28,      0
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
     -29,    -97,   -188,   -194,    186,   1444,   3886,   6365,   3881,   1441, 
     184,   -199,   -195,    -97,      1
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,    -29,    -97,   -188,   -194,    186,   1444,   3881,   6365,   3881, 
    1441,    184,   -199,   -195,    -96
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,    -29,    -97,   -188,   -194,    186,   1441,   3881,   6367, 
    3889,   1457,    199,   -215,   -315
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,    -29,    -97,   -189,   -195,    184,   1441,   3889, 
    6395,   3939,   1508,    146,   -609
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,    -29,    -98,   -192,   -199,    184,   1457, 
    3939,   6492,   4038,   1405,   -615
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0,    -30,   -101,   -195,   -199,    199, 
    1508,   4038,   6592,   3935,    640
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0,      0,    -28,    -97,   -195,   -215, 
     146,   1405,   3934,   6701,   4741
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0,      0,      0,      1,    -96,   -314, 
    -608,   -614,    638,   4732,  12637
}
};
int tik15_frame_cnt = 25;
int tik15_denominator = 16384;
int tik15_i_low[25] = {
       0,      0,      0,      0,      0,      0,      0,      0,      1,      2, 
       3,      4,      5,      6,      7,      8,      9,     10,     11,     12, 
      13,     14,     15,     16,     17
};
int tik15_i_high[25] = {
       7,      8,      9,     10,     11,     12,     13,     14,     15,     16, 
      17,     18,     19,     20,     21,     22,     23,     24,     24,     24, 
      24,     24,     24,     24,     24
};
int tik15_imatrix[25][25] = {
{

    8689,   5487,   2481,    647,   -153,   -341,   -273,   -151,      0,      0, 
       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

    5485,   5691,   3650,   1681,    459,    -84,   -220,   -180,   -101,      0, 
       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

    2481,   3652,   4885,   3460,   1749,    580,      9,   -170,   -161,    -99, 
       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

     647,   1682,   3461,   4952,   3582,   1842,    630,     27,   -167,   -165, 
    -103,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

    -153,    459,   1749,   3582,   5045,   3632,   1860,    632,     23,   -171, 
    -168,   -105,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

    -341,    -85,    580,   1842,   3632,   5063,   3634,   1856,    626,     20, 
    -173,   -168,   -105,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

    -273,   -220,      9,    630,   1860,   3634,   5059,   3630,   1847,    625, 
      20,   -173,   -168,   -105,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

    -151,   -180,   -170,     27,    632,   1856,   3630,   5058,   3616,   1848, 
     625,     20,   -172,   -168,   -105,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

       0,   -102,   -162,   -168,     23,    628,   1854,   3628,   5037,   3618, 
    1849,    625,     20,   -172,   -168,   -105,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

       0,      0,    -99,   -165,   -172,     20,    627,   1853,   3616,   5041, 
    3620,   1849,    625,     20,   -172,   -168,   -105,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

       0,      0,      0,   -104,   -168,   -173,     20,    627,   1847,   3618, 
    5044,   3620,   1849,    625,     20,   -173,   -168,   -105,      0,      0, 
       0,      0,      0,      0,      0
},
{

       0,      0,      0,      0,   -105,   -168,   -173,     20,    625,   1848, 
    3620,   5046,   3620,   1849,    626,     20,   -173,   -168,   -105,      0, 
       0,      0,      0,      0,      0
},
{

       0,      0,      0,      0,      0,   -105,   -168,   -173,     20,    625, 
    1849,   3620,   5046,   3620,   1850,    626,     20,   -173,   -168,   -105, 
       0,      0,      0,      0,      0
},
{

       0,      0,      0,      0,      0,      0,   -105,   -168,   -172,     20, 
     625,   1849,   3620,   5048,   3622,   1851,    626,     20,   -173,   -169, 
    -106,      0,      0,      0,      0
},
{

       0,      0,      0,      0,      0,      0,      0,   -105,   -168,   -172, 
      20,    625,   1849,   3620,   5045,   3624,   1851,    627,     20,   -174, 
    -171,   -108,      0,      0,      0
},
{

       0,      0,      0,      0,      0,      0,      0,      0,   -104,   -168, 
    -172,     20,    625,   1849,   3622,   5049,   3624,   1854,    626,     18, 
    -177,   -174,   -109,      0,      0
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,   -104, 
    -168,   -172,     20,    625,   1850,   3624,   5051,   3629,   1853,    625, 
      16,   -179,   -174,   -104,      0
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
    -104,   -168,   -172,     20,    626,   1851,   3624,   5056,   3629,   1855, 
     628,     19,   -178,   -182,   -128
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,   -104,   -168,   -173,     20,    625,   1851,   3629,   5061,   3638, 
    1869,    644,     25,   -217,   -319
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,   -105,   -168,   -174,     18,    624,   1855,   3637,   5079, 
    3673,   1910,    658,    -71,   -557
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,   -106,   -170,   -176,     16,    628,   1869,   3673, 
    5142,   3745,   1934,    491,   -663
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,   -108,   -173,   -179,     19,    644,   1911, 
    3746,   5231,   3773,   1736,   -212
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0,   -109,   -174,   -178,     25,    658, 
    1934,   3773,   5235,   3713,   1519
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0,      0,   -104,   -181,   -217,    -71, 
     490,   1734,   3709,   5698,   5343
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0,      0,      0,   -128,   -317,   -554, 
    -660,   -211,   1511,   5320,  11401
}
};
int tik20_frame_cnt = 25;
int tik20_denominator = 16384;
int tik20_i_low[25] = {
       0,      0,      0,      0,      0,      0,      0,      0,      1,      2, 
       3,      4,      5,      6,      7,      8,      9,     10,     11,     12, 
      13,     14,     15,     16,     17
};
int tik20_i_high[25] = {
       7,      8,      9,     10,     11,     12,     13,     14,     15,     16, 
      17,     18,     19,     20,     21,     22,     23,     24,     24,     24, 
      24,     24,     24,     24,     24
};
int tik20_imatrix[25][25] = {
{

    7604,   5343,   2910,   1159,    179,   -228,   -306,   -243,      0,      0, 
       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

    5321,   5205,   3591,   1929,    751,    101,   -164,   -210,   -162,      0, 
       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

    2899,   3591,   4224,   3181,   1850,    815,    197,    -85,   -160,   -139, 
       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

    1155,   1931,   3183,   4140,   3244,   1946,    894,    245,    -61,   -153, 
    -140,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

     179,    752,   1852,   3245,   4238,   3323,   1995,    917,    250,    -62, 
    -156,   -143,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

    -227,    101,    816,   1947,   3323,   4286,   3346,   2002,    908,    247, 
     -66,   -159,   -145,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

    -305,   -164,    197,    894,   1995,   3346,   4292,   3345,   1981,    904, 
     245,    -67,   -159,   -145,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

    -242,   -211,    -86,    246,    917,   2002,   3345,   4291,   3312,   1976, 
     903,    244,    -67,   -159,   -145,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

       0,   -164,   -162,    -62,    253,    916,   1999,   3342,   4250,   3307, 
    1977,    903,    244,    -67,   -159,   -145,      0,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

       0,      0,   -141,   -154,    -63,    249,    913,   1997,   3311,   4245, 
    3309,   1978,    904,    245,    -67,   -159,   -145,      0,      0,      0, 
       0,      0,      0,      0,      0
},
{

       0,      0,      0,   -141,   -158,    -66,    247,    912,   1978,   3307, 
    4247,   3311,   1979,    904,    245,    -67,   -159,   -146,      0,      0, 
       0,      0,      0,      0,      0
},
{

       0,      0,      0,      0,   -145,   -160,    -67,    246,    903,   1976, 
    3310,   4250,   3312,   1980,    904,    245,    -67,   -161,   -146,      0, 
       0,      0,      0,      0,      0
},
{

       0,      0,      0,      0,      0,   -146,   -161,    -68,    244,    902, 
    1978,   3311,   4250,   3313,   1980,    905,    245,    -68,   -161,   -147, 
       0,      0,      0,      0,      0
},
{

       0,      0,      0,      0,      0,      0,   -146,   -161,    -67,    244, 
     903,   1979,   3312,   4250,   3314,   1981,    905,    246,    -69,   -163, 
    -149,      0,      0,      0,      0
},
{

       0,      0,      0,      0,      0,      0,      0,   -146,   -159,    -67, 
     244,    904,   1979,   3313,   4255,   3315,   1981,    911,    245,    -71, 
    -165,   -152,      0,      0,      0
},
{

       0,      0,      0,      0,      0,      0,      0,      0,   -144,   -159, 
     -67,    244,    904,   1980,   3314,   4253,   3315,   1996,    911,    244, 
     -72,   -166,   -151,      0,      0
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,   -144, 
    -159,    -67,    244,    904,   1980,   3315,   4254,   3342,   1998,    913, 
     248,    -69,   -168,   -163,      0
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
    -144,   -159,    -67,    244,    904,   1981,   3316,   4291,   3348,   2008, 
     926,    258,    -74,   -213,   -282
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,   -145,   -160,    -68,    243,    904,   1982,   3348,   4302,   3371, 
    2037,    951,    247,   -175,   -479
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,   -146,   -162,    -70,    242,    906,   2008,   3371,   4342, 
    3421,   2078,    932,     77,   -625
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,   -148,   -164,    -71,    246,    927,   2038,   3421, 
    4407,   3474,   2054,    711,   -510
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,   -150,   -165,    -68,    258,    951,   2079, 
    3475,   4451,   3454,   1876,    234
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0,   -150,   -166,    -74,    247,    932, 
    2054,   3454,   4462,   3543,   2100
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0,      0,   -161,   -213,   -174,     77, 
     710,   1872,   3537,   5217,   5525
},
{

       0,      0,      0,      0,      0,      0,      0,      0,      0,      0, 
       0,      0,      0,      0,      0,      0,      0,   -281,   -477,   -622, 
    -508,    233,   2091,   5511,  10421
}
};
