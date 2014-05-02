
#include <iostream>
#include <glut.h>
#include <math.h>
#include <windows.h>                   

#define double float

double CL[4] = {(1 + sqrt(3.)) / (4 * sqrt(2.)),
                (3 + sqrt(3.)) / (4 * sqrt(2.)),
                (3 - sqrt(3.)) / (4 * sqrt(2.)),
                (1 - sqrt(3.)) / (4 * sqrt(2.))};
double CH[4] = {(1 - sqrt(3.)) / (4 * sqrt(2.)),
               -(3 - sqrt(3.)) / (4 * sqrt(2.)),
                (3 + sqrt(3.)) / (4 * sqrt(2.)),
               -(1 + sqrt(3.)) / (4 * sqrt(2.))}; 
typedef struct
{
  int W, H;
  double *Pix;
} PIC;

/* Glibal image and waveletted image */
PIC WorkPic, WavePic;

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

double* hpf_coeffs( double *Cl )
{
  int N = sizeof(Cl)/sizeof(*Cl);
  double* CHZ = (double *)malloc(N * sizeof(double));


  for (int k = 0; k < N; k++)
  {
    double j = pow(-1.0, k);
    CHZ[k] = j * Cl[N - k - 1];  
  }     

  return CHZ;
}

double* pconv( double *arr, int num, double *Cl, double *Ch)
{
  double* out = (double *)malloc(sizeof(double) * num);

  int M = 4;
  double sL, sH;

  for (int k = 0; k < num; k += 2)

  {
    sL = 0;
    sH = 0;
    for (int i = 0; i < 4; i++)
    {
      sL += arr[(k + i) % M] * Cl[i];
      sH += arr[(k + i) % M] * Ch[i];
    }
    out[k] = sL;
    out[k + 1] = sH;
  }
  return out;
}

PIC dwt2( PIC *P, double *Cl, double* Ch )
{
  PIC outpic, resoutpic;
  double *p;

  Create(&outpic, P->W, P->H);
  Create(&resoutpic, P->W, P->H);
  /* handle rows */

  for (int i = 0; i < P->H; i++)
  {
    p = pconv(&P->Pix[i * P->W], P->W, Cl, Ch);
    for (int j = 0; j < P->W; j++)
      outpic.Pix[i * P->W + j] = p[j]; //* 0.5;
    free(p);
  }

  /* handle columns*/

  for (int i = 0; /*0 && */i < P->W; i++)
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

  /* 0..h/2, 0..w/2 */
  for (int i = 0, k = 0; i < outpic.H / 2, k < outpic.H; i++, k += 2)
  {
    for (int j = 0, l = 0; j < outpic.W / 2, l < outpic.W; j++, l += 2)
    {
      resoutpic.Pix[i * outpic.W + j] = outpic.Pix[k * outpic.W + l];         
    }
  }

  /* h/2..h, 0..w/2 */
  for (int i = outpic.H / 2, k = 1; i < outpic.H, k < outpic.H; i++, k += 2)
  {
    for (int j = 0, l = 0; j < outpic.W / 2, l < outpic.W; j++, l += 2)
    {
      resoutpic.Pix[i * outpic.W + j] = outpic.Pix[k * outpic.W + l];         
    }
  }

  /* 0..h/2, w/2..w */
  for (int i = 0, k = 0; i < outpic.H / 2, k < outpic.H; i++, k += 2)
  {
    for (int j = outpic.W / 2, l = 1; j < outpic.W, l < outpic.W; j++, l += 2)
    {
      resoutpic.Pix[i * outpic.W + j] = outpic.Pix[k * outpic.W + l];         
    }
  }

  /* h/2..h, w/2..w */
  for (int i = outpic.H / 2, k = 1; i < outpic.H, k < outpic.H; i++, k += 2)
  {
    for (int j = outpic.W / 2, l = 1; j < outpic.W, l < outpic.W; j++, l += 2)
    {
      resoutpic.Pix[i * outpic.W + j] = outpic.Pix[k * outpic.W + l];         
    }
  }

  return resoutpic;
}

void Display( void )
{
  glClear(GL_COLOR_BUFFER_BIT);
  glPixelZoom(1, -1);
  glRasterPos2s(-1, 1);
  glDrawPixels(WorkPic.W, WorkPic.H, GL_LUMINANCE, GL_FLOAT, WorkPic.Pix);         //GL_DOUBLE
  glRasterPos2s(0, 1);
  glDrawPixels(WorkPic.W, WorkPic.H, GL_LUMINANCE, GL_FLOAT, WavePic.Pix);         //GL_DOUBLE

  glFinish();
  glutSwapBuffers();
}
void main( void )
{
  Load(&WorkPic, "4.bmp");
  WavePic = dwt2(&WorkPic, CL, CH);

  glutInitDisplayMode(GLUT_RGB);

  glutInitWindowPosition(0, 0);
  glutInitWindowSize(1000, 1000);
  glutCreateWindow("sdf");
  glutDisplayFunc(Display);
  glutMainLoop();
}
//void main( void )
//{
//  FILE *F;
//  BITMAPFILEHEADER fh;
//  BITMAPINFOHEADER ih;
//  RGBQUAD pal[1000];
//  unsigned char *pic, *row;
//  int bpl, c, r, g, b, x, y;
//
//  if ((F = fopen("3.bmp", "rb")) == NULL)
//  {
//    printf("wrong");
//    return;
//  }
//  if (fread(&fh, 1, sizeof(fh), F) != sizeof(fh))
//    printf("wrong header");
//
//  if (fread(&ih, 1, sizeof(ih), F) != sizeof(ih))
//    printf("wrong info header");  
//  
//  if (ih.biClrUsed == 0)
//  {
//    if (ih.biBitCount == 1)
//      c = 2;
//    else if (ih.biBitCount == 4)
//      c = 16;
//    else if (ih.biBitCount == 8)
//      c = 256;
//    else
//      c = 0;
//  }
//  else
//    c = ih.biClrUsed;
//
//  /* moved after headers, read the palette*/
//  fseek(F, fh.bfOffBits, SEEK_SET);
//  fread(pal, sizeof(RGBQUAD), c, F);
//
//  /* allocate memory now for picture, after palette */
//  pic = (unsigned char *)malloc(ih.biWidth * ih.biHeight);  
//  if (pic == NULL)
//    printf("memory problem\n");
//
//
//  if (ih.biBitCount == 1)
//    bpl = (ih.biWidth + 7) / 8;
//  else if (ih.biBitCount == 4)
//    bpl = (ih.biWidth + 1) / 2;
//  else if (ih.biBitCount == 8)
//    bpl = ih.biWidth;
//  else if (ih.biBitCount == 24)
//    bpl = ih.biWidth * 3;
//  else if (ih.biBitCount == 32)
//    bpl = ih.biWidth * 4;
//  else
//  {
//    fclose(F);
//    free(pic);
//  }
//
//  /* 4toby divided on 4 */
//  bpl = (bpl + 3) / 4 * 4;
//
//  if ((row = (unsigned char *)malloc(bpl)) == NULL)
//  {
//    printf("allocating error");
//    fclose(F);
//    free(pic);
//    return;
//  }
//
//  for (y = ih.biHeight - 1; y >= 0; y--)
//  {
//    fread(row, bpl, 1, F);
//    for (x = 0; x < ih.biWidth; x++)
//    {
//      if (ih.biBitCount == 1)
//      {
//        c = (row[x / 8] >> (7 - (x % 8))) & 1;
//        b = pal[c].rgbBlue;
//        g = pal[c].rgbGreen;
//        r = pal[c].rgbRed;
//      }
//      else if (ih.biBitCount == 4)
//      {
//        c = (row[x / 2] >> 4 * (1 - (x % 2))) & 0xF;
//        b = pal[c].rgbBlue;
//        g = pal[c].rgbGreen;
//        r = pal[c].rgbRed;        
//      }
//      else if (ih.biBitCount == 8)
//      {
//        c = row[x];
//        b = pal[c].rgbBlue;
//        g = pal[c].rgbGreen;
//        r = pal[c].rgbRed;
//      }
//      else if (ih.biBitCount == 24)
//      {
//        b = row[x * 3 + 0];
//        g = row[x * 3 + 1];
//        r = row[x * 3 + 2];
//      }
//      else if (ih.biBitCount == 32)
//      {
//        b = row[x * 4 + 0];
//        g = row[x * 4 + 1];
//        r = row[x * 4 + 2];
//      }
//      c = (r * 30 + g * 59 + b * 11) / 100;
//      pic[y * ih.biWidth + x] = c;                        
//    }
//  } 
//
//  free(row);
//  fclose(F);
//
//  for (int i = 0; i < 4; i++)
//    printf("%.16f\n", CH[i]);
//
//  
//
//  W = ih.biWidth;
//  H = ih.biHeight;
//  Img = pic;
//
//  double C[2] = {0.5, 0.5};
//  double CC[2] = {0.5, -0.5};
//  double corn[2] = {0.5, 0.5};
//  double arrray[4] = {1, 2, 3, 4};
//  dwt2(arrray, C, CC, 4, 1);
//
//
//
//  glutInitDisplayMode(GLUT_RGB);
//  glutInitWindowPosition(0, 0);
//  glutInitWindowSize(1000, 1000);
//  glutCreateWindow("sdf");
//  glutDisplayFunc(Display);
//  glutMainLoop();
//
//  free(pic);
//} 

