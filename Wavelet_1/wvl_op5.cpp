#include <iostream>
#include <glut.h>
#include <math.h>
#include <windows.h>                   
#include <stdlib.h>

#define double float

//double CL[2] = {0.5, 0.5};
//double CH[2] = {0.5, -0.5};

float MIN, MAX;

double CL[4] = {(1 + sqrt(3.)) / (4 * sqrt(2.)),    // c0
                (3 + sqrt(3.)) / (4 * sqrt(2.)),    // c1
                (3 - sqrt(3.)) / (4 * sqrt(2.)),    // c2
                (1 - sqrt(3.)) / (4 * sqrt(2.))};   // c3
double CH[4] = {(1 - sqrt(3.)) / (4 * sqrt(2.)),    // c3
               -(3 - sqrt(3.)) / (4 * sqrt(2.)),    // -c2
                (3 + sqrt(3.)) / (4 * sqrt(2.)),    // c1
               -(1 + sqrt(3.)) / (4 * sqrt(2.))};   // -c0
double OCL[4] = {(3 - sqrt(3.)) / (4 * sqrt(2.)),   // c2
                 (3 + sqrt(3.)) / (4 * sqrt(2.)),   // c1
                 (1 + sqrt(3.)) / (4 * sqrt(2.)),   // c0
                 (1 - sqrt(3.)) / (4 * sqrt(2.))};  // c3                
double OCH[4] = {(1 - sqrt(3.)) / (4 * sqrt(2.)),   // c3                
                -(1 + sqrt(3.)) / (4 * sqrt(2.)),   // -c0
                 (3 + sqrt(3.)) / (4 * sqrt(2.)),   // c1
                -(3 - sqrt(3.)) / (4 * sqrt(2.))};  // -c2
typedef struct
{
  int W, H;
  double *Pix;
} PIC;

/* Global image and waveletted image */
PIC WorkPic, WavePic;

/* Min value */
float Min( float *data, int N )
{
  float minval = data[0];

  for (int i = 0; i < N - 1; i++)
  {
    if (minval > data[i + 1])
      minval = data[i + 1];
  }

  return minval;
}

float Max( float *data, int N )
{
  float maxval;

  maxval = data[0];
  for (int i = 0; i < N - 1; i++)
  {
    if (maxval < data[i+1])
      maxval = data[i+1];
  }

  return maxval;
}

/* Creating picture: allocate memory, initialize W, H */ 
int Create( PIC *P, int W, int H )
{
  P->W = P->H = 0;
  P->Pix = (double *)malloc(H * W * sizeof(double));
  if (P->Pix == NULL)
    return 0;

  P->H = H;
  P->W = W;
  return 1;
}

/* Load picture from file using filename */
int Load( PIC *P, char *FileName )
{
  FILE *F;
  BITMAPFILEHEADER fh;
  BITMAPINFOHEADER ih;
  RGBQUAD pal[1000];
  unsigned char *row;
  int bpl, c, r, g, b, x, y;

  P->W = P->H = 0;
  P->Pix = NULL;
  if ((F = fopen(FileName, "rb")) == NULL)
  {
    printf("wrong");
    return 0;
  }
  if (fread(&fh, 1, sizeof(fh), F) != sizeof(fh))
  {
    fclose(F);
    printf("wrong header");
    return 0;
  }
  if (fread(&ih, 1, sizeof(ih), F) != sizeof(ih))
  {
    printf("wrong info header");  
    fclose(F);
    return 0;
  }
  
  if (ih.biClrUsed == 0)
  {
    if (ih.biBitCount == 1)
      c = 2;
    else if (ih.biBitCount == 4)
      c = 16;
    else if (ih.biBitCount == 8)
      c = 256;
    else
      c = 0;
  }
  else
    c = ih.biClrUsed;

  /* moved after headers, read the palette*/
  fseek(F, fh.bfOffBits, SEEK_SET);
  fread(pal, sizeof(RGBQUAD), c, F);

  if (ih.biBitCount == 1)
    bpl = (ih.biWidth + 7) / 8;
  else if (ih.biBitCount == 4)
    bpl = (ih.biWidth + 1) / 2;
  else if (ih.biBitCount == 8)
    bpl = ih.biWidth;
  else if (ih.biBitCount == 24)
    bpl = ih.biWidth * 3;
  else if (ih.biBitCount == 32)
    bpl = ih.biWidth * 4;
  else
  {
    fclose(F);
    return 0;
  }

  /* 4toby divided on 4 */
  bpl = (bpl + 3) / 4 * 4;

  if ((row = (unsigned char *)malloc(bpl)) == NULL)
  {
    printf("allocating error");
    fclose(F);
    return 0;
  }

  if (!Create(P, ih.biWidth, ih.biHeight))
  {
    free(row);
    fclose(F);
    return 0;
  }

  for (y = ih.biHeight - 1; y >= 0; y--)
  {
    fread(row, bpl, 1, F);
    for (x = 0; x < ih.biWidth; x++)
    {
      if (ih.biBitCount == 1)
      {
        c = (row[x / 8] >> (7 - (x % 8))) & 1;
        b = pal[c].rgbBlue;
        g = pal[c].rgbGreen;
        r = pal[c].rgbRed;
      }
      else if (ih.biBitCount == 4)
      {
        c = (row[x / 2] >> 4 * (1 - (x % 2))) & 0xF;
        b = pal[c].rgbBlue;
        g = pal[c].rgbGreen;
        r = pal[c].rgbRed;        
      }
      else if (ih.biBitCount == 8)
      {
        c = row[x];
        b = pal[c].rgbBlue;
        g = pal[c].rgbGreen;
        r = pal[c].rgbRed;
      }
      else if (ih.biBitCount == 24)
      {
        b = row[x * 3 + 0];
        g = row[x * 3 + 1];
        r = row[x * 3 + 2];
      }
      else if (ih.biBitCount == 32)
      {
        b = row[x * 4 + 0];
        g = row[x * 4 + 1];
        r = row[x * 4 + 2];
      }
      c = (r * 30 + g * 59 + b * 11) / 100;
      P->Pix[y * ih.biWidth + x] = c / 255.0;                        
    }
  } 

  free(row);
  fclose(F);
  return 1;
}
         
/* Coefficient to high-freq-filter
 * NOT USED */
double * hpf_coeffs( double *Cl )
{
  int N = sizeof(Cl)/sizeof(*Cl);
  double *CHZ = (double *)malloc(N * sizeof(double));


  for (int k = 0; k < N; k++)
  {
    double j = pow(-1.0, k);
    CHZ[k] = j * Cl[N - k - 1];  
  }     

  return CHZ;
}

/* Pair convolution. Dimensional transform */
double *pconv( double *arr, int num, double *Cl, double *Ch)
{
  double *out = (double *)malloc(sizeof(double) * num);

  int M = 2;
  double sL, sH;

  for (int k = 0; k < num; k += 2)
  {
    sL = 0;
    sH = 0;
    for (int i = 0; i < 4; i++)
    {
      sL += arr[(k + i) % num] * Cl[i];
      sH += arr[(k + i) % num] * Ch[i];
    }
    out[k] = sL;
    out[k + 1] = sH;
  }
  return out;
}

double *ipconv( double *arr, int num, double *Cl, double *Ch)
{
  double *out = (double *)malloc(sizeof(double) * num);

  int M = 2;
  double sL, sH;
  bool flag = false;

  for (int k = num - 2; k < num; k += 2)
  {
    sL = 0;
    sH = 0;
    for (int i = 0; i < 4; i++)
    {
      sL += arr[(k + i) % num] * Cl[i];
      sH += arr[(k + i) % num] * Ch[i];
    }
    out[k] = sL;
    out[k + 1] = sH;
    if (k == (num - 2) && flag == false)
    {
      k = 0;
      flag = true;
    }
  }
  return out;
}


/* 2-dimensional transform */
PIC dwt2( PIC *P, double *Cl, double* Ch )
{
  PIC outpic, resoutpic;
  double *p;
  double min, max;

  Create(&outpic, P->W, P->H);
  Create(&resoutpic, P->W, P->H);
  
  /* handle columns*/      
  for (int i = 0; i < P->W; i++)
  {
    double *v = (double *)malloc(sizeof(double) * P->H);
    for (int j = 0; j < P->H; j++)
      v[j] = P->Pix[j * P->W + i];
    p = pconv(v, P->H, Cl, Ch);
    for (int j = 0; j < P->H; j++)
      outpic.Pix[j * P->W + i] = p[j];
    free(p);
    free(v);
  }

  /* handle rows */
  for (int i = 0; i < P->H; i += 2)
  {
    p = pconv(&P->Pix[i * P->W], P->W, Cl, Ch);
    for (int j = 0; j < P->W; j++)
      outpic.Pix[i * P->W + j] = p[j]; //* 0.5;  
    free(p);
  }                      
  
  
  MIN = Min(outpic.Pix, outpic.H * outpic.W);
  MAX = Max(outpic.Pix, outpic.H * outpic.W);
  /*
  printf("Min\n");
  printf("%.16f\n", Min(outpic.Pix, outpic.H * outpic.W));
  printf("Max\n");

  printf("%.16f\n", Max(outpic.Pix, outpic.H * outpic.W));
  */

  /* 0..h/2, 0..w/2 */
  for (int i = 0, k = 0; i < outpic.H / 2 && k < outpic.H; i++, k += 2)
    for (int j = 0, l = 0; j < outpic.W / 2 &&  l < outpic.W; j++, l += 2)
      resoutpic.Pix[i * outpic.W + j] = outpic.Pix[k * outpic.W  + l];         

  /*
  for (int i = 0; i < outpic.H; i++)
    for (int j = 0; j < outpic.W; j++)
      resoutpic.Pix[i * outpic.W + j] = outpic.Pix[i * outpic.W + j];
  */
  /* h/2..h, 0..w/2 */
  for (int i = outpic.H / 2, k = 1; i < outpic.H && k < outpic.H; i++, k += 2)
    for (int j = 0, l = 0; j < outpic.W / 2 && l < outpic.W; j++, l += 2)
      resoutpic.Pix[i * outpic.W + j] = outpic.Pix[k * outpic.W + l];         
  /* 0..h/2, w/2..w */
  for (int i = 0, k = 0; i < outpic.H / 2 && k < outpic.H; i++, k += 2)
    for (int j = outpic.W / 2, l = 1; j < outpic.W && l < outpic.W; j++, l += 2)
      resoutpic.Pix[i * outpic.W + j] = outpic.Pix[k * outpic.W + l];         

  /* h/2..h, w/2..w */
  for (int i = outpic.H / 2, k = 1; i < outpic.H && k < outpic.H; i++, k += 2)
    for (int j = outpic.W / 2, l = 1; j < outpic.W && l < outpic.W; j++, l += 2)
      resoutpic.Pix[i * outpic.W + j] = outpic.Pix[k * outpic.W + l];         


  for (int i = 0; i < resoutpic.W * resoutpic.H; i++)
    resoutpic.Pix[i] = (resoutpic.Pix[i] - MIN) / (MAX - MIN);

  return resoutpic;
}

PIC idwt2( PIC *P, double *iCl, double *iCh )
{
  PIC resoutpic, outpic;   
  float min, max;
  double *p;

  Create(&outpic, P->W, P->H);
  Create(&resoutpic, P->W, P->H);
  //min = Min(P->Pix, P->H * P->W);
  //max = Max(P->Pix, P->H * P->W);
  /*
  printf("Min\n");
  printf("%.16f\n", Min(outpic.Pix, outpic.H * outpic.W));
  printf("Max\n");

  printf("%.16f\n", Max(outpic.Pix, outpic.H * outpic.W));
  */

  /* 0..h/2, 0..w/2 */
  for (int i = 0, k = 0; k < P->H; i++, k += 2)
    for (int j = 0, l = 0; l < P->W; j++, l += 2)
      outpic.Pix[k * outpic.W + l] = P->Pix[i * outpic.W  + j];         

  /* h/2..h, 0..w/2 */
  for (int i = P->H / 2, k = 1; i < P->H && k < P->H; i++, k += 2)
    for (int j = 0, l = 0; j < P->W / 2 && l < P->W; j++, l += 2)
      outpic.Pix[k * outpic.W + l] = P->Pix[i * outpic.W  + j];         

  /* 0..h/2, w/2..w */
  for (int i = 0, k = 0; i < P->H / 2 && k < P->H; i++, k += 2)
    for (int j = P->W / 2, l = 1; j < P->W && l < P->W; j++, l += 2)
      outpic.Pix[k * outpic.W + l] = P->Pix[i * outpic.W  + j];               

  /* h/2..h, w/2..w */
  for (int i = P->H / 2, k = 1; i < P->H && k < P->H; i++, k += 2)
    for (int j = P->W / 2, l = 1; j < P->W && l < P->W; j++, l += 2)
      outpic.Pix[k * outpic.W + l] = P->Pix[i * outpic.W  + j];               


  for (int i = 0; i < outpic.W * outpic.H; i++)
    outpic.Pix[i] = (outpic.Pix[i] + MIN) * (MAX - MIN);

  /* handle columns*/      
  for (int i = 0; i < outpic.W; i++)
  {
    double *v = (double *)malloc(sizeof(double) * outpic.H);
    for (int j = 0; j < outpic.H; j++)
      v[j] = outpic.Pix[j * outpic.W + i];
    p = ipconv(v, resoutpic.H, iCl, iCh);
    for (int j = 0; j < outpic.H; j++)
      resoutpic.Pix[j * outpic.W + i] = p[j];
    free(p);
    free(v);
  }               

  /* handle rows */
  for (int i = 0; i < outpic.H; i += 2)
  {
    p = ipconv(&outpic.Pix[i * outpic.W], outpic.W, iCl, iCh);
    for (int j = 0; j < outpic.W; j++)
      resoutpic.Pix[i * outpic.W + j] = p[j]; //* 0.5;  
    free(p);
  }  

  min = Min(resoutpic.Pix, resoutpic.H * resoutpic.W);
  max = Max(resoutpic.Pix, resoutpic.H * resoutpic.W);
/*
  for (int i = 0; i < resoutpic.W * resoutpic.H; i++)
    resoutpic.Pix[i] = (resoutpic.Pix[i] - min) / (max - min);*/


  return resoutpic;
}

void Display( void )
{
  glClearColor(1, 0, 0, 1);
  glClear(GL_COLOR_BUFFER_BIT);
  glPixelZoom(1, -1);
  glRasterPos2s(-1, 1);
  glDrawPixels(WorkPic.W, WorkPic.H, GL_LUMINANCE, GL_FLOAT, WorkPic.Pix);         //GL_DOUBLE
  glRasterPos2s(0, 1);
  glDrawPixels(WorkPic.W, WorkPic.H, GL_LUMINANCE, GL_FLOAT, WavePic.Pix);         //GL_DOUBLE

  glFinish();
  glutSwapBuffers();
}

void Keyboard( unsigned char Key, int X, int Y )
{
  if (Key == 27)
    exit(0);
  else if (Key == 'f')
    glutFullScreen();
}

void main( void )
{
  //for (int i = 0; i < 4; i++)
    //printf("%.16f ", CL[i]);//*sqrt(2.0));
  Load(&WorkPic, "3.bmp");
  WavePic = dwt2(&WorkPic, CL, CH);
  WorkPic = idwt2(&WavePic, OCL, OCH);
  glutInitDisplayMode(GLUT_RGB);

  glutInitWindowPosition(0, 0);
  glutInitWindowSize(1000, 1000);
  glutCreateWindow("sdf");

  glutDisplayFunc(Display);
  glutKeyboardFunc(Keyboard);

  glutMainLoop();
}
