#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <math.h>
#include <stdlib.h>

#include <windows.h>

#include <glut.h>

#define MAIN_WINDOW_TITLE "sdf:Wavelet Sample"

int WinWidth, WinHeight;
char ImageFileName[1000] = "D:\\4.bmp";

void Reshape( int Width, int Height )
{
  glViewport(0, 0, WinWidth = Width, WinHeight = Height);
} /* End of 'Reshape' function */


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
PIC WorkPic, WavePic, ResPic;

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
  W = (W + 1) / 2 * 2;
  H = (H + 1) / 2 * 2;
  P->Pix = (double *)malloc(H * W * sizeof(double));
  if (P->Pix == NULL)
    return 0;
  memset(P->Pix, 0, H * W * sizeof(double));
  P->H = H;
  P->W = W;
  return 1;
}

/* Creating picture: allocate memory, initialize W, H */
void FreePic( PIC *P )
{
  if (P->Pix != NULL)
    free(P->Pix);
  P->Pix = NULL;
  P->W = P->H = 0;
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

/* Pair convolution. Dimensional transform */
double *pconv( double *arr, int num, double *Cl, double *Ch)
{
  double *out = (double *)malloc(sizeof(double) * num);
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

  for (int k = num - 2, j = 0; k < num && j < num; k += 2, j += 2)
  {
    sL = 0;
    sH = 0;
    for (int i = 0; i < 4; i++)
    {
      sL += arr[(k + i) % num] * Cl[i];
      sH += arr[(k + i) % num] * Ch[i];
    }
    out[j] = sL;
    out[j + 1] = sH;
    if (k == (num - 2) && flag == false)
    {
      k = -2;
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

  Create(&outpic, P->W, P->H);
  Create(&resoutpic, P->W, P->H);

  for (int i = 0; i < P->W * P->H; i++)
    outpic.Pix[i] = 0.5;


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
  for (int i = 0; i < P->H; i++)
  {
    p = pconv(&outpic.Pix[i * P->W], P->W, Cl, Ch);
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

  /* 0..h/2, 0..w/2 */
  for (int i = 0, k = 0; k < P->H; i++, k += 2)
    for (int j = 0, l = 0; l < P->W; j++, l += 2)
      outpic.Pix[k * outpic.W + l] = P->Pix[i * outpic.W  + j];

  /* h/2..h, 0..w/2 */
  for (int i = P->H / 2, k = 1; i < P->H; i++, k += 2)
    for (int j = 0, l = 0; j < P->W / 2; j++, l += 2)
      outpic.Pix[k * outpic.W + l] = P->Pix[i * outpic.W  + j];

  /* 0..h/2, w/2..w */
  for (int i = 0, k = 0; i < P->H / 2; i++, k += 2)
    for (int j = P->W / 2, l = 1; j < P->W; j++, l += 2)
      outpic.Pix[k * outpic.W + l] = P->Pix[i * outpic.W  + j];

  /* h/2..h, w/2..w */
  for (int i = P->H / 2, k = 1; i < P->H; i++, k += 2)
    for (int j = P->W / 2, l = 1; j < P->W; j++, l += 2)
      outpic.Pix[k * outpic.W + l] = P->Pix[i * outpic.W  + j];

  for (int i = 0; i < outpic.W * outpic.H; i++)
    outpic.Pix[i] = outpic.Pix[i] * (MAX - MIN) + MIN;

  /* handle rows */
  for (int i = 0; i < outpic.H; i++)
  {
    p = ipconv(&outpic.Pix[i * outpic.W], outpic.W, iCl, iCh);
    for (int j = 0; j < outpic.W; j++)
      resoutpic.Pix[i * outpic.W + j] = p[j];
    free(p);
  }

  /* handle columns */
  for (int i = 0; i < outpic.W; i++)
  {
    double *v = (double *)malloc(sizeof(double) * outpic.H);
    for (int j = 0; j < outpic.H; j++)
      v[j] = resoutpic.Pix[j * outpic.W + i];
    p = ipconv(v, resoutpic.H, iCl, iCh);
    for (int j = 0; j < outpic.H; j++)
      resoutpic.Pix[j * outpic.W + i] = p[j];
    free(p);
    free(v);
  }

  min = Min(resoutpic.Pix, resoutpic.H * resoutpic.W);
  max = Max(resoutpic.Pix, resoutpic.H * resoutpic.W);

  for (int i = 0; i < resoutpic.W * resoutpic.H; i++)
    resoutpic.Pix[i] = (resoutpic.Pix[i] - min) / (max - min);

  return resoutpic;
}

void Display( void )
{
  static char Txt[] = "'f' - full screen | 'o' - open image | 'esc' - exit";
  static char Buf[1000];

  glClearColor(1, 1, 1, 1);
  glClear(GL_COLOR_BUFFER_BIT);

  glLoadIdentity();
  //glRasterPos2s(-0.98, -0.98);
  glColor3d(1, 1, 1);
  glListBase(1111);
  glRasterPos2d(-0.999, -0.99);
  glCallLists(sprintf(Buf, "%s", Txt), GL_UNSIGNED_BYTE, Buf);
  glRasterPos2d(-0.999, -0.88);
  glCallLists(sprintf(Buf, "file name: \"%s\"", ImageFileName),
GL_UNSIGNED_BYTE, Buf);

  glPixelZoom(1, -1);
  glRasterPos2d(-1, 1);
  glDrawPixels(WorkPic.W, WorkPic.H, GL_LUMINANCE, GL_FLOAT,
WorkPic.Pix);         // GL_DOUBLE
  glRasterPos2d(-0.333, 1);
  glDrawPixels(WavePic.W, WavePic.H, GL_LUMINANCE, GL_FLOAT,
WavePic.Pix);         // GL_DOUBLE
  glRasterPos2d(0.333, 1);
  glDrawPixels(ResPic.W, ResPic.H, GL_LUMINANCE, GL_FLOAT,
ResPic.Pix);         // GL_DOUBLE

  glFinish();
  glutSwapBuffers();
  glutPostRedisplay();
}

void Idle( void )
{
  glutPostRedisplay();
}

float entropy( float *data, int n )
{
  float *arr = new float[n];
  float *p = new float[n];
  int kol = 0;
  float a;
  int count;
  bool flag;
  float entropy = 0;


  for (int i = 0; i < n; i++)
  {
    flag = true;
    a = data[i];
    count = 1;

    for (int k = 0 ; k < kol; k++)
    {
      if (a == arr[k])
      {
        flag = false;
        break;
      }
    }

    if (flag == false)
      continue;

    for (int j = i + 1; j < n; j++)
    {
      if (a == data[j])
        count++;
    }
    arr[kol] = a;
    p[kol++] = count;
  }
  return kol;
  for (int i = 0; i < kol; i++)
    entropy -= (p[i] / n) * log(p[i] / n) / log(2.);

  delete [] arr;
  delete [] p;

  return entropy;   
}

void Keyboard( unsigned char Key, int X, int Y )
{
  static int IsFullScreen = 0, SaveW, SaveH;

  if (Key == 27)
    exit(0);
  else if (Key == 'f')
  {
    if (IsFullScreen)
      /* Восстанавливаем размеры окна */
      glutReshapeWindow(SaveW, SaveH);
    else
    {
      SaveW = WinWidth;
      SaveH = WinHeight;
      /* Переходим в полноэкранный режим */
      glutFullScreen();
    }
    IsFullScreen = !IsFullScreen;
  }
  else if (Key == 'o')
  {
    HWND hWnd = FindWindow(NULL, (TCHAR *)MAIN_WINDOW_TITLE);
    static OPENFILENAME ofn = {0};

    ofn.lStructSize = sizeof(OPENFILENAME);
    ofn.Flags = OFN_FILEMUSTEXIST;
    ofn.hwndOwner = hWnd;
    ofn.lpstrFile = ImageFileName;
    ofn.nMaxFile = sizeof(ImageFileName);
    ofn.lpstrTitle = "Open Image File";
    ofn.lpstrFilter = 
      "All supported image formats\0*.BMP\0"
      "Microsoft Bitmap\0*.BMP"
      "All files (*.*)\0*.*\0";
    if (GetOpenFileName(&ofn))
    {
      FreePic(&WorkPic);
      FreePic(&WavePic);
      FreePic(&ResPic);

      /* Loading image */
      //if (Load(&WorkPic, ImageFileName))
      //{
      //  WavePic = dwt2(&WorkPic, CL, CH);
      //  for (int i = 0; i < WorkPic.H / 2 ; i++)
      //    for (int j =  0; j < WorkPic.W / 2; j++)
      //    {
      //      double x = WavePic.Pix[(i + WorkPic.H / 2) * WorkPic.W + j + WorkPic.W / 2];
      //      double y = WavePic.Pix[i * WorkPic.W + j];
      //      double z = WavePic.Pix[i * WorkPic.W + j + WorkPic.W / 2];
      //      double p = WavePic.Pix[(i + WorkPic.H / 2) * WorkPic.W + j]; 
      //      x = (int)(x * 4) / 4.0;
      //      y = (int)(y * 32) / 32.0;
      //      z = (int)(z * 32) / 32.0;
      //      p = (int)(z * 16) / 16.0;

      //      WavePic.Pix[(i + WorkPic.H / 2) * WorkPic.W + j + WorkPic.W / 2] = x;
      //      //WavePic.Pix[i * WorkPic.W + j] = y;
      //      //WavePic.Pix[i * WorkPic.W + j + WorkPic.W / 2] = z;
      //      //WavePic.Pix[(i + WorkPic.H / 2) * WorkPic.W + j] = p;
      //      //WavePic.Pix[i * WorkPic.W + j] = 0;
      //    }
     
        
       // std :: cout << entropy(WorkPic.Pix, WorkPic.W * WorkPic.H)<<"\n";
       /* ResPic = idwt2(&WavePic, OCL, OCH);
        std :: cout << entropy(ResPic.Pix, WorkPic.W * WorkPic.H)<<"\n";
      } */   
      Load(&WorkPic, ImageFileName);
      std :: cout << entropy(WorkPic.Pix, WorkPic.W * WorkPic.H)<<"\n";
      WavePic = dwt2(&WorkPic, CL, CH);
      ResPic = idwt2(&WavePic, OCL, OCH); 
    }
  }
}



void main( void )
{
  HDC hDC;

  static float rurk[1000000];

  glutInitDisplayMode(GLUT_RGB);
  glutInitWindowPosition(0, 0);
  glutInitWindowSize(1000, 1000);
  glutCreateWindow(MAIN_WINDOW_TITLE);

  glutDisplayFunc(Display);
  glutKeyboardFunc(Keyboard);
  glutIdleFunc(Idle);
  glutReshapeFunc(Reshape);

  hDC = GetDC(NULL);
  SelectObject(hDC,
    CreateFont(30, 0, 0, 0, FW_BOLD, 0, 0, 0, RUSSIAN_CHARSET,
      OUT_TT_PRECIS, CLIP_DEFAULT_PRECIS, PROOF_QUALITY, FF_MODERN | FIXED_PITCH, ""));
  wglUseFontBitmaps(hDC, 0, 256, 1111);
  ReleaseDC(NULL, hDC);

  glutMainLoop();
}
