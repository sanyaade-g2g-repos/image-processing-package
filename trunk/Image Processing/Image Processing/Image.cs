using System;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.IO;
using System.Drawing.Imaging;
using System.Collections;
using System.Windows.Forms;
using EngMATLib;
using MatlabMultiplication;

namespace Image_Task1
{
    public class Image
    {
        #region Member Variables
        //Private Members
        private string path;
        private string type;
        private string comment;
        private Size size;
        private string range;
        private int[] redhistogram;
        private int[] greenhistogram;
        private int[] bluehistogram;
        private double[] newredhistogram;
        private double[] newgreenhistogram;
        private double[] newbluehistogram;
        private int redmax = 0;
        private int greenmax = 0;
        private int bluemax = 0;
        int MinR, MaxR, MinG, MaxG, MinB, MaxB;
        private Bitmap histogram;
        private bool IsModified = false;
        
        double[,] RedReal;
        double[,] GreenReal;
        double[,] BlueReal;
        double[,] RedImg;
        double[,] GreenImg; 
        double[,] BlueImg ;
        double[,] MagR;
        double[,] MagG;
        double[,] MagB;

        //Public Members
        public Bitmap BitmapImage;
        public Bitmap Fourier;

        //For Unsafe Code
        public struct PixelData
        {
            public byte blue;
            public byte green;
            public byte red;
        }
        private BitmapData bmData;
        private IntPtr Scan0;
        private int stride;
        #endregion

        #region Constructors
        public Image()
        {
        }

        public Image(Bitmap b)
        {
            BitmapImage = b;
        }

        public Image(string _path)
        {
            path = _path;
            IsModified = false;
        }

        public Image(string _path, string _type, string _comment, Size _size,string _range)
        {
            path = _path;
            type = _type;
            comment = _comment;
            size = _size;
            IsModified = false;
        }
        #endregion

        public Bitmap OpenImage()
        {
            if (Path.GetExtension(path).ToLower() == ".ppm")
            {
                StreamReader sr = new StreamReader(path);
                type = sr.ReadLine();
                comment = sr.ReadLine();
                string _size = sr.ReadLine();
                string[] tokens = _size.Split(' ');
                size.Width = int.Parse(tokens[0]);
                size.Height = int.Parse(tokens[1]);

                range = sr.ReadLine();

                BitmapImage = new Bitmap(size.Width, size.Height);
                BitmapData BitmapImageData1 = BitmapImage.LockBits(new Rectangle(0, 0, BitmapImage.Width, BitmapImage.Height),
                    ImageLockMode.WriteOnly, PixelFormat.Format32bppArgb);

                int i = 0, j = 0, k = 0, m = 0;
                int[] rgb = new int[3];
                int strideOffset = BitmapImageData1.Stride - (BitmapImageData1.Width * 4);

                unsafe
                {
                    byte* imgPtr1 = (byte*)BitmapImageData1.Scan0;
                    while (!sr.EndOfStream)
                    {
                        k = 0;
                        tokens = sr.ReadLine().Split(' ');
                        while (k < tokens.Length)
                        {
                            if (tokens[k] != "")
                                rgb[m++] = int.Parse(tokens[k]);

                            k++;
                            if (m == 3)
                            {
                                if (i < size.Width && j < size.Height)
                                {
                                    imgPtr1[0] = (byte)rgb[2]; // B
                                    imgPtr1[1] = (byte)rgb[1]; // G
                                    imgPtr1[2] = (byte)rgb[0]; // R
                                    imgPtr1[3] = 255;           // A 

                                    i++;
                                    imgPtr1 += 4;
                                }
                                if (i >= size.Width)
                                {
                                    j++;
                                    i = 0;
                                    imgPtr1 += strideOffset;
                                }
                                if (j >= size.Height) break;
                                m = 0;
                            }
                        }
                    }
                }

                BitmapImage.UnlockBits(BitmapImageData1);
                sr.Close();
            }
            else
            {
                BitmapImage = new Bitmap(path);
            }
            this.InitializeHistogram();
            return BitmapImage;
        }

        #region Image Geometry

        public Image ImageResize(int new_width, int new_height)
        {
            int width = BitmapImage.Width;
            int height = BitmapImage.Height;
            double wratio = (double)new_width * 1.0 / (double)width;
            double hratio = (double)new_height * 1.0 / (double)height;
            int i, j;
            Image new_img = new Image("temp.ppm", "p3", "After Resizing", new Size(new_width, new_height), "100");
            new_img.BitmapImage = new Bitmap(new_width, new_height);

            //unsafe
            {
                #region Resizing Code
                for (j = 0; j < new_height; j++)
                {
                    double y = (double)(j / hratio);
                    int ypixel = (int)y;
                    double dy = y - (double)ypixel;

                    for (i = 0; i < new_width; i++)
                    {
                        double x = (double)(i / wratio);
                        int xpixel = (int)x;
                        double dx = x - (double)xpixel;

                        if (xpixel == width - 1 || ypixel == height - 1)
                            continue;

                        Color c1 = BitmapImage.GetPixel(xpixel, ypixel);
                        Color c2 = BitmapImage.GetPixel(xpixel, ypixel + 1);
                        Color c3 = BitmapImage.GetPixel(xpixel + 1, ypixel);
                        Color c4 = BitmapImage.GetPixel(xpixel + 1, ypixel + 1);

                        double xInterpolated1_R = (double)((double)c1.R * (1 - dy) + (double)c2.R * (dy));
                        double xInterpolated1_G = (double)((double)c1.G * (1 - dy) + (double)c2.G * (dy));
                        double xInterpolated1_B = (double)((double)c1.B * (1 - dy) + (double)c2.B * (dy));

                        double xInterpolated2_R = (double)((double)c3.R * (1 - dy) + (double)c4.R * (dy));
                        double xInterpolated2_G = (double)((double)c3.G * (1 - dy) + (double)c4.G * (dy));
                        double xInterpolated2_B = (double)((double)c3.B * (1 - dy) + (double)c4.B * (dy));

                        double pixelvalue_R = xInterpolated1_R * (1 - dx) + xInterpolated2_R * (dx);
                        double pixelvalue_G = xInterpolated1_G * (1 - dx) + xInterpolated2_G * (dx);
                        double pixelvalue_B = xInterpolated1_B * (1 - dx) + xInterpolated2_B * (dx);

                        //Slow Code
                        new_img.BitmapImage.SetPixel(i, j, Color.FromArgb((int)pixelvalue_R, (int)pixelvalue_G, (int)pixelvalue_B));

                        //Fast Code
                        //imgPtr1[0] = (int)pixelvalue;
                        //imgPtr1++;
                        /*imgPtr1[0] = ;     // B
                        imgPtr1[1] = g;     // G
                        imgPtr1[2] = r;     // R
                        imgPtr1[3] = 255;   // A 
                        imgPtr1 += 4;*/

                    }
                    //imgPtr1 += strideOffset;
                }
                #endregion
            }
            this.BitmapImage = new_img.BitmapImage;
            this.size = new_img.size;
            return new_img;
        }

        #endregion

        #region Pixels Operations

        public void InitializeHistogram()
        {
            int i, j, pix;
            redhistogram = new int[256];
            greenhistogram = new int[256];
            bluehistogram = new int[256];

            for (i = 0; i < 256; i++)
            {
                redhistogram[i] = 0;
                greenhistogram[i] = 0;
                bluehistogram[i] = 0;
            }

            for (i = 0; i < BitmapImage.Width; i++)
            {
                for (j = 0; j < BitmapImage.Height; j++)
                {
                    Color c = BitmapImage.GetPixel(i, j);
                    pix = c.R;
                    redhistogram[pix]++;
                    if (redhistogram[pix] > redmax)
                    {
                        redmax = redhistogram[pix];
                    }

                    pix = c.G;
                    greenhistogram[pix]++;
                    if (greenhistogram[pix] > greenmax)
                    {
                        greenmax = greenhistogram[pix];
                    }

                    pix = c.B;
                    bluehistogram[pix]++;
                    if (bluehistogram[pix] > bluemax)
                    {
                        bluemax = bluehistogram[pix];
                    }
                }
            }
        }

        public Bitmap Histogram(Color c)
        {
            histogram = new Bitmap(400, 300);

            int[] dummy = new int[256];
            int i, j;
            int max = redmax;

            redhistogram.CopyTo(dummy, 0);
            if (c == Color.Green)
            {
                greenhistogram.CopyTo(dummy, 0);
                max = greenmax;
            }
            else if (c == Color.Blue)
            {
                bluehistogram.CopyTo(dummy, 0);
                max = bluemax;
            }

            for (i = 255; i >= 0; i--)
            {
                for (j = 299; j >= 0; j--)
                {
                    if (dummy[i] > 0)
                    {
                        histogram.SetPixel(i, j, c);

                        //To fit PictureBox
                        dummy[i] -= (max / 250);
                    }
                }
            }
            return this.histogram;
        }

        public Image EditContrast(int value)
        {
            Image Edited = new Image();
            Edited.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);
            double[,] red = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            double[,] blue = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            double[,] green = new double[this.BitmapImage.Width, this.BitmapImage.Height];

            Color OldMin = GetMinVal(BitmapImage);
            Color OldMax = GetMaxVal(BitmapImage);

            for (int i = 0; i < BitmapImage.Width; i++)
                for (int j = 0; j < BitmapImage.Height; j++)
                {
                    Color OldPixel = BitmapImage.GetPixel(i, j);
                    red[i, j] = (int)   ((((((float)OldPixel.R / 255.0) - 0.5) * value/100.0) + 0.5) * 255);
                    green[i, j] = (int) ((((((float)OldPixel.G / 255.0) - 0.5) * value/100.0) + 0.5) * 255);
                    blue[i, j] = (int)  ((((((float)OldPixel.B / 255.0) - 0.5) * value/100.0) + 0.5) * 255);

                    Color NewPixel = CutOff((int)red[i, j], (int)green[i, j], (int)blue[i, j]);
                    Edited.BitmapImage.SetPixel(i, j, NewPixel);
                }
            return Edited;
        }

        public Image EditBrightness(int Degree)
        {
            Image Edited = new Image();
            Edited.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);

            for (int i = 0; i < BitmapImage.Width; i++)
                for (int j = 0; j < BitmapImage.Height; j++)
                {
                    Color OldPixel = BitmapImage.GetPixel(i, j);

                    int R = OldPixel.R + Degree;
                    int G = OldPixel.G + Degree;
                    int B = OldPixel.B + Degree;

                    Color NewPixel = CutOff(R, G, B);
                    Edited.BitmapImage.SetPixel(i, j, NewPixel);
                }
            return Edited;
        }

        public Image ImageAlgebra(Image First, Image Second, string Operation)
        {
            // Resizing Second  to First
            Second.ImageResize(First.BitmapImage.Width, First.BitmapImage.Height);
            Image NewImage = new Image(new Bitmap(First.BitmapImage.Width, First.BitmapImage.Height));
            //System.Windows.Forms.MessageBox.Show(First.BitmapImage.GetPixel(0,0).ToString());
            //System.Windows.Forms.MessageBox.Show(Second.BitmapImage.GetPixel(0, 0).ToString());
            //int OldMinA = AllOperations(First.BitmapImage.GetPixel(0, 0).A, Second.BitmapImage.GetPixel(0, 0).A, Operation);
            //int OldMaxA = AllOperations(First.BitmapImage.GetPixel(0, 0).A, Second.BitmapImage.GetPixel(0, 0).A, Operation);
            int OldMinR = AllOperations(First.BitmapImage.GetPixel(0, 0).R, Second.BitmapImage.GetPixel(0, 0).R, Operation);
            int OldMaxR = AllOperations(First.BitmapImage.GetPixel(0, 0).R, Second.BitmapImage.GetPixel(0, 0).R, Operation);
            int OldMinG = AllOperations(First.BitmapImage.GetPixel(0, 0).G, Second.BitmapImage.GetPixel(0, 0).G, Operation);
            int OldMaxG = AllOperations(First.BitmapImage.GetPixel(0, 0).G, Second.BitmapImage.GetPixel(0, 0).G, Operation);
            int OldMinB = AllOperations(First.BitmapImage.GetPixel(0, 0).B, Second.BitmapImage.GetPixel(0, 0).B, Operation);
            int OldMaxB = AllOperations(First.BitmapImage.GetPixel(0, 0).B, Second.BitmapImage.GetPixel(0, 0).B, Operation);
            
            // Getting The Min & Max Value in The New Bitmap
            for (int i = 0; i < NewImage.BitmapImage.Width; i++)
                for (int j = 0; j < NewImage.BitmapImage.Height; j++)
                {
                    Color OldPixel1 = First.BitmapImage.GetPixel(i, j);
                    Color OldPixel2 = Second.BitmapImage.GetPixel(i, j);

                    //int A = AllOperations(OldPixel1.A, OldPixel2.A, Operation);
                    int R = AllOperations(OldPixel1.R, OldPixel2.R, Operation);
                    int G = AllOperations(OldPixel1.G, OldPixel2.G, Operation);
                    int B = AllOperations(OldPixel1.B, OldPixel2.B, Operation);

                    //if (A > OldMaxA) OldMaxA = A;
                    //else if (A < OldMinA) OldMinA = A;
                    if (R > OldMaxR) OldMaxR = R;
                    else if (R < OldMinR) OldMinR = R;
                    if (G > OldMaxG) OldMaxG = G;
                    else if (G < OldMinG) OldMinG = G;
                    if (B > OldMaxB) OldMaxB = B;
                    else if (B < OldMinB) OldMinB = B;
                }

            // Normalizing The New Bitmap
            for (int i = 0; i < NewImage.BitmapImage.Width; i++)
                for (int j = 0; j < NewImage.BitmapImage.Height; j++)
                {
                    Color OldPixel1 = First.BitmapImage.GetPixel(i, j);
                    Color OldPixel2 = Second.BitmapImage.GetPixel(i, j);

                    //int A = AllOperations(OldPixel1.A, OldPixel2.A, Operation);
                    int R = System.Math.Abs(AllOperations(OldPixel1.R, OldPixel2.R, Operation));
                    int G = System.Math.Abs(AllOperations(OldPixel1.G, OldPixel2.G, Operation));
                    int B = System.Math.Abs(AllOperations(OldPixel1.B, OldPixel2.B, Operation));

                    //if (R < 0) R *= -1; if (G < 0) G *= -1; if (B < 0) B *= -1;
                    //if (G < 0) System.Windows.Forms.MessageBox.Show(G.ToString());

                    //A = Normalize(A, OldMinA, OldMaxA, 0, 255);
                    R = Normalize(R, OldMinR, OldMaxR, 0, 255);
                    G = Normalize(G, OldMinG, OldMaxG, 0, 255);
                    B = Normalize(B, OldMinB, OldMaxB, 0, 255);

                    try{ NewImage.BitmapImage.SetPixel(i, j, Color.FromArgb(R, G, B)); }
                    catch {  }
                }
            this.BitmapImage = NewImage.BitmapImage;
            return NewImage;
        }

        public Image EditQuantization(int Degree)
        {
            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);

            int NoOfBits = (int)System.Math.Log(Degree, 2);
            int Begin = 7;
            int ToBeAnded = 0;
            for (int i = 0; i < NoOfBits; i++)
            {
                ToBeAnded += (int)System.Math.Pow(2, Begin);
                Begin--;
            }
           
            for (int i = 0; i < BitmapImage.Width; i++)
                for (int j = 0; j < BitmapImage.Height; j++)
                {
                    Color OldPixel = BitmapImage.GetPixel(i, j);

                    byte R = (byte)OldPixel.R;
                    byte G = (byte)OldPixel.G;
                    byte B = (byte)OldPixel.B;

                    int NR = (int)(R & ToBeAnded);
                    int NG = (int)(G & ToBeAnded);
                    int NB = (int)(B & ToBeAnded);

                    NewImage.BitmapImage.SetPixel(i, j, Color.FromArgb(NR, NG, NB));
                }

            return NewImage;
        }
                
        public Image Negative()
        {
            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);

            for (int i = 0; i < NewImage.BitmapImage.Width; i++)
                for (int j = 0; j < NewImage.BitmapImage.Height; j++)
                {
                    NewImage.BitmapImage.SetPixel(i, j, Color.FromArgb(255 - this.BitmapImage.GetPixel(i, j).R, 255 - this.BitmapImage.GetPixel(i, j).G, 255 - this.BitmapImage.GetPixel(i, j).B));
                }

            return NewImage;
        }

        public Image Power(double pow)
        {
            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);
            double[,] BufferR = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            double[,] BufferG = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            double[,] BufferB = new double[this.BitmapImage.Width, this.BitmapImage.Height];

            for (int i = 0; i < NewImage.BitmapImage.Width; i++)
                for (int j = 0; j < NewImage.BitmapImage.Height; j++)
                {
                    BufferR[i, j] = Math.Pow((double)this.BitmapImage.GetPixel(i, j).R + 1, pow);
                    BufferG[i, j] = Math.Log((double)this.BitmapImage.GetPixel(i, j).G + 1, pow);
                    BufferB[i, j] = Math.Log((double)this.BitmapImage.GetPixel(i, j).B + 1, pow);
                }

            double MinR = GetMinVal(BufferR); double MaxR = GetMaxVal(BufferR);
            double MinG = GetMinVal(BufferG); double MaxG = GetMaxVal(BufferG);
            double MinB = GetMinVal(BufferB); double MaxB = GetMaxVal(BufferB);

            NewImage.BitmapImage = Normalize(BufferR, BufferB, BufferG, 0, 255);

            return NewImage;
        }

        public Image Log()
        {
            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);
            double[,] BufferR = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            double[,] BufferG = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            double[,] BufferB = new double[this.BitmapImage.Width, this.BitmapImage.Height];

            for (int i = 0; i < NewImage.BitmapImage.Width; i++)
                for (int j = 0; j < NewImage.BitmapImage.Height; j++)
                {
                    BufferR[i, j] = Math.Log((double)this.BitmapImage.GetPixel(i, j).R + 1);
                    BufferG[i, j] = Math.Log((double)this.BitmapImage.GetPixel(i, j).G + 1);
                    BufferB[i, j] = Math.Log((double)this.BitmapImage.GetPixel(i, j).B + 1);
                }

            double MinR = GetMinVal(BufferR); double MaxR = GetMaxVal(BufferR);
            double MinG = GetMinVal(BufferG); double MaxG = GetMaxVal(BufferG);
            double MinB = GetMinVal(BufferB); double MaxB = GetMaxVal(BufferB);

            NewImage.BitmapImage = Normalize(BufferR, BufferB, BufferG, 0, 255);

            return NewImage;
        }

        #endregion

        #region Illuminaion

        public Image SingleScaleRetinex(double Const, Image OldImage)
        {
            double n = (3.7 * Const) - 0.5;
            int MaskSize = (int)((2 * n) + 1.0);
            if (MaskSize % 2 == 0) MaskSize++;

            double[,] Mask = new double[(int)MaskSize, (int)MaskSize];
            double[,] BlurR = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            double[,] BlurG = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            double[,] BlurB = new double[this.BitmapImage.Width, this.BitmapImage.Height];

            Image NewImage = new Image();

            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);

            for (int a = -MaskSize / 2; a <= MaskSize / 2; a++)
            {
                for (int b = -MaskSize / 2; b <= MaskSize / 2; b++)
                {
                    double X = (2.0 * Math.PI) * Math.Pow(Const, 2);
                    double Y = 1.0 / X;
                    double powe = -1 * (Math.Pow(a, 2) + Math.Pow(b, 2));
                    double exp = powe / (2.0 * Math.Pow(Const, 2));
                    double Z = Y * (Math.Exp(exp));
                    Mask[(int)(a + MaskSize / 2), (int)(b + MaskSize / 2)] = Z;
                }
            }

            //Convulution
            int border = (int)(MaskSize / 2.0);
            for (int i = 0; i < BitmapImage.Width; i++)
            {
                for (int j = 0; j < BitmapImage.Height; j++)
                {
                    double NewR = 0, NewG = 0, NewB = 0;

                    for (int m = -1 * border; m <= border; m++)
                    {
                        for (int k = -1 * border; k <= border; k++)
                        {
                            if ((m + i) >= 0 && (m + i) < BitmapImage.Width && (k + j) >= 0 && (k + j) < BitmapImage.Height)
                            {
                                Color c = BitmapImage.GetPixel(m + i, k + j);
                                NewR += Mask[m + border, k + border] * (double)c.R;
                                NewG += Mask[m + border, k + border] * (double)c.G;
                                NewB += Mask[m + border, k + border] * (double)c.B;
                            }
                        }
                    }
                    BlurR[i, j] = NewR;
                    BlurG[i, j] = NewG;
                    BlurB[i, j] = NewB;
                }
            }

            for (int i = 0; i < BitmapImage.Width; i++)
                for (int j = 0; j < BitmapImage.Height; j++)
                {
                    Color c = OldImage.BitmapImage.GetPixel(i, j);

                    BlurR[i, j] = Math.Log((double)(c.R + 1.0) / (double)(BlurR[i, j] + 1.0));
                    BlurG[i, j] = Math.Log((double)(c.G + 1.0) / (double)(BlurG[i, j] + 1.0));
                    BlurB[i, j] = Math.Log((double)(c.B + 1.0) / (double)(BlurB[i, j] + 1.0));
                }

            Normalizedouble(BlurR, BlurG, BlurB);
            for (int t = 0; t < BitmapImage.Width; t++)
            {
                for (int j = 0; j < BitmapImage.Height; j++)
                {
                    NewImage.BitmapImage.SetPixel(t, j, Color.FromArgb((int)BlurR[t, j], (int)BlurG[t, j], (int)BlurB[t, j]));
                }
            }
            return NewImage;

        }

        public Image GammaIntensityCorrection(int C, double G)
        {
            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);
            double[,] NewR = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            double[,] NewB = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            double[,] NewG = new double[this.BitmapImage.Width, this.BitmapImage.Height];

            for (int i = 0; i < this.BitmapImage.Width; i++)
            {
                for (int j = 0; j < this.BitmapImage.Height; j++)
                {
                    NewR[i,j] = C * (Math.Pow((double)(BitmapImage.GetPixel(i, j).R), 1.0 / (double)(G + 0.0001)));
                    NewG[i, j] = C * (Math.Pow((double)(BitmapImage.GetPixel(i, j).G), 1.0 / (double)(G + 0.0001)));
                    NewB[i, j] = C * (Math.Pow((double)(BitmapImage.GetPixel(i, j).B), 1.0 / (double)(G + 0.0001)));
                }
            }
            NewImage.BitmapImage = Normalize(NewR, NewG, NewB, 0, 255);
            return NewImage;
        }

        public Image HistogramEqualization()
        {
            InitializeHistogram();

            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(BitmapImage.Width, BitmapImage.Height);

            newredhistogram = new double[256];
            newgreenhistogram = new double[256];
            newbluehistogram = new double[256];

            redhistogram.CopyTo(newredhistogram, 0);
            bluehistogram.CopyTo(newbluehistogram, 0);
            greenhistogram.CopyTo(newgreenhistogram, 0);

            double RedMax;
            double GreenMax;
            double BlueMax;

            for (int i = 1; i <= 255; i++)
            {
                newredhistogram[i] += newredhistogram[i - 1];
                newgreenhistogram[i] += newgreenhistogram[i - 1];
                newbluehistogram[i] += newbluehistogram[i - 1];
            }

            RedMax = newredhistogram[255];
            GreenMax = newgreenhistogram[255];
            BlueMax = newbluehistogram[255];

            for (int i = 0; i <= 255; i++)
            {
                newredhistogram[i] /= RedMax;
                newgreenhistogram[i] /= GreenMax;
                newbluehistogram[i] /= BlueMax;
            }

            for (int i = 0; i <= 255; i++)
            {
                newredhistogram[i] *= 255;
                newgreenhistogram[i] *= 255;
                newbluehistogram[i] *= 255;
            }

            for (int i = 0; i < BitmapImage.Width; i++)
            {
                for (int j = 0; j < BitmapImage.Height; j++)
                {
                    double NewR = 0; double NewG = 0; double NewB = 0;

                    int newpix;
                    Color c = BitmapImage.GetPixel(i, j);
                    newpix = c.R;

                    for (int k = 0; k <= 255; k++)
                    {
                        if (newpix == k)
                        {
                            NewR = newredhistogram[k];
                            break;
                        }
                    }

                    newpix = c.G;
                    for (int k = 0; k <= 255; k++)
                    {
                        if (newpix == k)
                        {
                            NewG = newgreenhistogram[k];
                            break;
                        }
                    }

                    newpix = c.B;
                    for (int k = 0; k <= 255; k++)
                    {
                        if (newpix == k)
                        {
                            NewB = newbluehistogram[k];
                            break;
                        }
                    }

                    try { NewImage.BitmapImage.SetPixel(i, j, CutOff((int)NewR, (int)NewG, (int)NewB)); }
                    catch { }
                }
               
            }
            NewImage.InitializeHistogram();
            return NewImage;
        }

        public Image HistogramMatching(Image Image2)
        {
            Image NewImage = new Image();

            Image2.ImageResize(this.BitmapImage.Width, this.BitmapImage.Height);
            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);

            double Rpixno, Gpixno, Bpixno;
            int indexr = 0;
            double minr1 = 0, minr2 = 0;
            int indexg = 0;
            double ming1 = 0, ming2 = 0;
            int indexb = 0;
            double minb1 = 0, minb2 = 0;

            double[] NewRedhistogram = new double[256];
            double[] NewGreenhistogram = new double[256];
            double[] NewBluehistogram = new double[256];

            double[] FirstRedhistogram = new double[256];
            double[] FirstGreenhistogram = new double[256];
            double[] FirstBluehistogram = new double[256];

            double[] SecRedhistogram = new double[256];
            double[] SecGreenhistogram = new double[256];
            double[] SecBluehistogram = new double[256];

            //this.InitializeHistogram();
            this.HistogramEqualization();
            newredhistogram.CopyTo(FirstRedhistogram, 0);
            newgreenhistogram.CopyTo(FirstGreenhistogram, 0);
            newbluehistogram.CopyTo(FirstBluehistogram, 0);

            //Image2.InitializeHistogram();
            Image2.HistogramEqualization();
            Image2.newredhistogram.CopyTo(SecRedhistogram, 0);
            Image2.newgreenhistogram.CopyTo(SecGreenhistogram, 0);
            Image2.newbluehistogram.CopyTo(SecBluehistogram, 0);

            for (int i = 0; i <= 255; i++)
            {
                Rpixno = FirstRedhistogram[i];
                Gpixno = FirstGreenhistogram[i];
                Bpixno = FirstBluehistogram[i];

                for (int j = 0; j <= 255; j++)
                {
                    if (Rpixno == SecRedhistogram[j])
                    {
                        NewRedhistogram[i] = j;
                        break;
                    }
                    else
                    {
                        for (int k = 0; k <= 255; k++)
                        {
                            if (k == 0)
                            {
                                minr1 = Math.Abs((double)FirstRedhistogram[i] - (double)SecRedhistogram[k]);
                                indexr = i;
                            }
                            else
                            {
                                minr2 = Math.Abs((double)FirstRedhistogram[i] - (double)SecRedhistogram[k]);

                                if (minr2 < minr1)
                                {
                                    minr1 = minr2;
                                    indexr = k;
                                }
                            }
                        }
                        NewRedhistogram[i] = indexr;

                    }
                }
                for (int j = 0; j <= 255; j++)
                {
                    if (Gpixno == SecGreenhistogram[j])
                    {
                        NewGreenhistogram[i] = j;
                        break;
                    }
                    else
                    {
                        for (int k = 0; k <= 255; k++)
                        {
                            if (k == 0)
                            {
                                ming1 = Math.Abs((double)FirstGreenhistogram[i] - (double)SecGreenhistogram[k]);
                                indexg = i;
                            }
                            else
                            {
                                ming2 = Math.Abs((double)FirstGreenhistogram[i] - (double)SecGreenhistogram[k]);

                                if (ming2 < ming1)
                                {
                                    ming1 = ming2;
                                    indexg = k;
                                }
                            }
                        }
                        NewGreenhistogram[i] = indexg;

                    }
                }
                for (int j = 0; j <= 255; j++)
                {
                    if (Bpixno == SecBluehistogram[j])
                    {
                        NewBluehistogram[i] = j;
                        break;
                    }
                    else
                    {
                        for (int k = 0; k <= 255; k++)
                        {
                            if (k == 0)
                            {
                                minb1 = Math.Abs((double)FirstBluehistogram[i] - (double)SecBluehistogram[k]);
                                indexb = i;
                            }
                            else
                            {
                                minb2 = Math.Abs((double)FirstBluehistogram[i] - (double)SecBluehistogram[k]);

                                if (minb2 < minb1)
                                {
                                    minb1 = minb2;
                                    indexb = k;
                                }
                            }
                        }
                        NewBluehistogram[i] = indexb;

                    }
                }
            } int NewR = 0, NewG = 0, NewB = 0;

            /**********************************Displaying Image******************************/
            for (int X = 0; X < NewImage.BitmapImage.Width; X++)
                for (int Y = 0; Y < NewImage.BitmapImage.Height; Y++)
                {
                    Color c = this.BitmapImage.GetPixel(X, Y);
                    for (int a = 0; a <= 255; a++)
                        if (c.R == a)
                        {
                            NewR = (int)NewRedhistogram[a];
                            break;
                        }
                    for (int a = 0; a <= 255; a++)
                        if (c.G == a)
                        {
                            NewG = (int)NewGreenhistogram[a];
                            break;
                        }
                    for (int a = 0; a <= 255; a++)
                        if (c.B == a)
                        {
                            NewB = (int)NewBluehistogram[a];
                            break;
                        }
                    { NewImage.BitmapImage.SetPixel(X, Y, Color.FromArgb(NewR, NewG, NewB)); }
                    //catch { }
                }
            /*********************************************************************************/
            return NewImage;
        }

        #endregion Illumination

        #region Filters
        public Image Convulution(double[,]Mask)
        {
            Image NewImage = new Image();
            int Border = (Mask.GetLength(0) - 1) / 2;
            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);

            bmData = this.BitmapImage.LockBits(new Rectangle(0, 0, this.BitmapImage.Width, this.BitmapImage.Height),
                ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
            stride = bmData.Stride;
            Scan0 = bmData.Scan0;

            unsafe
            {
                double[,] BufferR = new double[this.BitmapImage.Width, this.BitmapImage.Height];
                double[,] BufferG = new double[this.BitmapImage.Width, this.BitmapImage.Height];
                double[,] BufferB = new double[this.BitmapImage.Width, this.BitmapImage.Height];

                for (int i = Border; i < this.BitmapImage.Width - Border; i++)
                    for (int j = Border; j < this.BitmapImage.Height - Border; j++)
                    {
                        PixelData* pPixel = (PixelData*)(Scan0.ToInt32() + j * stride + i * sizeof(PixelData));

                        BufferR[i, j] = pPixel->red;
                        BufferG[i, j] = pPixel->green;
                        BufferB[i, j] = pPixel->blue;
                    }

                double MinR = GetMinVal(BufferR); double MaxR = GetMaxVal(BufferR);
                double MinG = GetMinVal(BufferG); double MaxG = GetMaxVal(BufferG);
                double MinB = GetMinVal(BufferB); double MaxB = GetMaxVal(BufferB);

                for (int i = Border; i < NewImage.BitmapImage.Height - Border; i++)
                    for (int j = Border; j < NewImage.BitmapImage.Width - Border; j++)
                    {
                        double NewR = 0; double NewG = 0; double NewB = 0;

                        for (int m = -1 * Border; m <= Border; m++)
                            for (int k = -1 * Border; k <= Border; k++)
                            {
                                PixelData* pPixel = (PixelData*)(Scan0.ToInt32() + (i + k) * stride + (j + m) * sizeof(PixelData));
                                NewR += (Mask[m + Border, k + Border] * pPixel->red);
                                NewG += (Mask[m + Border, k + Border] * pPixel->green);
                                NewB += (Mask[m + Border, k + Border] * pPixel->blue);
                            }

                        BufferR[j, i] = NewR;
                        BufferG[j, i] = NewG;
                        BufferB[j, i] = NewB;
                    }
                NewImage.BitmapImage = Normalize(BufferR, BufferG, BufferB, 0, 255);
            }
            this.BitmapImage.UnlockBits(bmData);
            return NewImage;
        }
        
        public Image Blur(int Length)
        {
            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);
            int Border = (Length - 1) / 2;

            for (int i = Border; i < NewImage.BitmapImage.Width - Border; i++)
                for (int j = Border; j < NewImage.BitmapImage.Height - Border; j++)
                {
                    int NewR = 0; int NewG = 0; int NewB = 0;

                    for (int m = -1 * Border; m <= Border; m++)
                        for (int k = -1 * Border; k <= Border; k++)
                        {
                            NewR += (BitmapImage.GetPixel(i + m, j + k).R);
                            NewG += (BitmapImage.GetPixel(i + m, j + k).G);
                            NewB += (BitmapImage.GetPixel(i + m, j + k).B);
                        }
                    NewR /= Length * Length; NewG /= Length * Length; NewB /= Length * Length;
                    NewImage.BitmapImage.SetPixel(i, j, Color.FromArgb(NewR, NewG, NewB));
                }

            return NewImage;            
        }

        public Image Detection(int[,] Mask)
        {
            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);

            Color OldMin = GetMinVal(BitmapImage);
            Color OldMax = GetMaxVal(BitmapImage);

            int MinR = 0; int MinG = 0; int MinB = 0;
            int MaxR = 0; int MaxG = 0; int MaxB = 0;

            for (int i = 1; i < NewImage.BitmapImage.Width - 1; i++)
                for (int j = 1; j < NewImage.BitmapImage.Height - 1; j++)
                {
                    int NewR = 0; int NewG = 0; int NewB = 0;

                    for (int m = -1; m <= 1; m++)
                        for (int k = -1; k <= 1; k++)
                        {
                            NewR += (Mask[m + 1, k + 1] * BitmapImage.GetPixel(i + m, j + k).R);
                            NewG += (Mask[m + 1, k + 1] * BitmapImage.GetPixel(i + m, j + k).G);
                            NewB += (Mask[m + 1, k + 1] * BitmapImage.GetPixel(i + m, j + k).B);
                        }

                    if (NewR < MinR)MinR = NewR;
                    else if (NewR > MaxR) MaxR = NewR;
                    if (NewG < MinG)MinG = NewG;
                    else if (NewG > MaxG) MaxG = NewG;
                    if (NewR < MinB)MinB = NewB;
                    else if (NewB > MaxB) MaxB = NewB;
                }
 
                for (int i = 1; i < NewImage.BitmapImage.Width - 1; i++)
                    for (int j = 1; j < NewImage.BitmapImage.Height - 1; j++)
                    {
                        int NewR = 0; int NewG = 0; int NewB = 0;
                        for (int m = -1; m <= 1; m++)
                            for (int k = -1; k <= 1; k++)
                            {
                                NewR += (Mask[m + 1, k + 1] * BitmapImage.GetPixel(i + m, j + k).R);
                                NewG += (Mask[m + 1, k + 1] * BitmapImage.GetPixel(i + m, j + k).G);
                                NewB += (Mask[m + 1, k + 1] * BitmapImage.GetPixel(i + m, j + k).B);
                            }

                        try { NewImage.BitmapImage.SetPixel(i, j, Color.FromArgb(Normalize(NewR, MinR, MaxR, 0, 255), Normalize(NewG, MinG, MaxG, 0, 255), Normalize(NewB, MinB, MaxB, 0, 255))); }
                        catch {  }
                    }

            return NewImage;
        }

        public Image Sharpening(int[,]Mask)
        {
            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);

            Color OldMin = GetMinVal(BitmapImage);
            Color OldMax = GetMaxVal(BitmapImage);

            for (int i = 1; i < NewImage.BitmapImage.Width - 1; i++)
                for (int j = 1; j < NewImage.BitmapImage.Height - 1; j++)
                {
                    int NewR = 0; int NewG = 0; int NewB = 0;
                    for (int m = -1; m <= 1; m++)
                        for (int k = -1; k <= 1; k++)
                        {
                            NewR += (Mask[m + 1, k + 1] * BitmapImage.GetPixel(i + m, j + k).R);
                            NewG += (Mask[m + 1, k + 1] * BitmapImage.GetPixel(i + m, j + k).G);
                            NewB += (Mask[m + 1, k + 1] * BitmapImage.GetPixel(i + m, j + k).B);
                        }

                    try { NewImage.BitmapImage.SetPixel(i, j, CutOff(NewR, NewG, NewB)); }
                    catch { }
                }

            return NewImage;
        }
        
        #endregion

        #region Frequency Domain

        public Image IdealLPF(double Do)
        {
            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);
            NewImage.Fourier = this.Fourier;

            NewImage.RedImg = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];
            NewImage.GreenImg = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];
            NewImage.BlueImg = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];

            NewImage.RedReal = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];
            NewImage.GreenReal = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];
            NewImage.BlueReal = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];

            for (int u = 0; u < this.RedImg.GetLength(0); u++)
                for (int v = 0; v < this.RedImg.GetLength(1); v++)
                {
                    double D = Math.Sqrt(Math.Pow((u - (this.RedImg.GetLength(0)*1.0 / 2.0)), 2) + Math.Pow((v - (this.RedImg.GetLength(1)*1.0 / 2.0)), 2));
                    if (D > Do)
                    {
                        NewImage.RedImg[u, v] = 0;
                        NewImage.GreenImg[u, v] = 0;
                        NewImage.BlueImg[u, v] = 0;

                        NewImage.RedReal[u, v] = 0;
                        NewImage.GreenReal[u, v] = 0;
                        NewImage.BlueReal[u, v] = 0;
                    }
                    else
                    {
                        NewImage.RedImg[u, v] = this.RedImg[u, v];
                        NewImage.GreenImg[u, v] = this.GreenImg[u, v];
                        NewImage.BlueImg[u, v] = this.BlueImg[u, v];

                        NewImage.RedReal[u, v] = this.RedReal[u, v];
                        NewImage.GreenReal[u, v] = this.GreenReal[u, v];
                        NewImage.BlueReal[u, v] = this.BlueReal[u, v];
                    }
                }

            return NewImage;
        }

        public Image ButterWorthLPF(double Do, double N)
        {
            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);
            NewImage.Fourier = this.Fourier;

            NewImage.RedImg = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];
            NewImage.GreenImg = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];
            NewImage.BlueImg = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];

            NewImage.RedReal = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];
            NewImage.GreenReal = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];
            NewImage.BlueReal = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];

            for (int u = 0; u < this.RedImg.GetLength(0); u++)
                for (int v = 0; v < this.RedImg.GetLength(1); v++)
                {
                    double D = Math.Sqrt(Math.Pow((u - (this.RedImg.GetLength(0) * 1.0 / 2.0)), 2) + Math.Pow((v - (this.RedImg.GetLength(1) * 1.0 / 2.0)), 2));

                    NewImage.RedImg[u, v] = (1.0 / (1.0 + Math.Pow(D / Do, 2.0 * N))) * this.RedImg[u, v];
                    NewImage.GreenImg[u, v] = (1.0 / (1.0 + Math.Pow(D / Do, 2.0 * N))) * this.GreenImg[u, v];
                    NewImage.BlueImg[u, v] = (1.0 / (1.0 + Math.Pow(D / Do, 2.0 * N))) * this.BlueImg[u, v];

                    NewImage.RedReal[u, v] = (1.0 / (1.0 + Math.Pow(D / Do, 2.0 * N))) * this.RedReal[u, v];
                    NewImage.GreenReal[u, v] = (1.0 / (1.0 + Math.Pow(D / Do, 2.0 * N))) * this.GreenReal[u, v];
                    NewImage.BlueReal[u, v] = (1.0 / (1.0 + Math.Pow(D / Do, 2.0 * N))) * this.BlueReal[u, v];
                }

            return NewImage;
        }

        public Image GaussianLPF(double Gamma)
        {
            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);
            NewImage.Fourier = this.Fourier;

            NewImage.RedImg = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];
            NewImage.GreenImg = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];
            NewImage.BlueImg = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];

            NewImage.RedReal = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];
            NewImage.GreenReal = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];
            NewImage.BlueReal = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];

            for (int u = 0; u < this.BitmapImage.Width; u++)
                for (int v = 0; v < this.BitmapImage.Height; v++)
                {
                    double D = Math.Sqrt(Math.Pow((u - (this.BitmapImage.Width * 1.0 / 2.0)), 2) + Math.Pow((v - (this.BitmapImage.Height * 1.0 / 2.0)), 2));

                    double NewR = Math.Exp(-1.0 * Math.Pow(D, 2.0) / (2.0 * Math.Pow(Gamma, 2))) * this.RedImg[u, v];
                    double NewG = Math.Exp(-1.0 * Math.Pow(D, 2.0) / (2.0 * Math.Pow(Gamma, 2))) * this.GreenImg[u, v];
                    double NewB = Math.Exp(-1.0 * Math.Pow(D, 2.0) / (2.0 * Math.Pow(Gamma, 2))) * this.BlueImg[u, v];

                    NewImage.RedImg[u, v] = NewR;
                    NewImage.GreenImg[u, v] = NewG;
                    NewImage.BlueImg[u, v] = NewB;
                }

            for (int u = 0; u < this.BitmapImage.Width; u++)
                for (int v = 0; v < this.BitmapImage.Height; v++)
                {
                    double D = Math.Sqrt(Math.Pow((u - (this.BitmapImage.Width * 1.0 / 2.0)), 2) + Math.Pow((v - (this.BitmapImage.Height * 1.0 / 2.0)), 2));

                    double NewR = Math.Exp(-1.0 * Math.Pow(D, 2.0) / (2.0 * Math.Pow(Gamma, 2))) * this.RedReal[u, v];
                    double NewG = Math.Exp(-1.0 * Math.Pow(D, 2.0) / (2.0 * Math.Pow(Gamma, 2))) * this.GreenReal[u, v];
                    double NewB = Math.Exp(-1.0 * Math.Pow(D, 2.0) / (2.0 * Math.Pow(Gamma, 2))) * this.BlueReal[u, v];

                    NewImage.RedReal[u, v] = NewR;
                    NewImage.GreenReal[u, v] = NewG;
                    NewImage.BlueReal[u, v] = NewB;
                }

            return NewImage;
        }

        public Image IdealHPF(double Do)
        {
            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);
            NewImage.Fourier = this.Fourier;

            NewImage.RedImg = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];
            NewImage.GreenImg = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];
            NewImage.BlueImg = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];

            NewImage.RedReal = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];
            NewImage.GreenReal = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];
            NewImage.BlueReal = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];

            for (int u = 0; u < this.RedImg.GetLength(0); u++)
                for (int v = 0; v < this.RedImg.GetLength(1); v++)
                {
                    double D = Math.Sqrt(Math.Pow((u - (this.RedImg.GetLength(0) / 2)), 2) + Math.Pow((v - (this.RedImg.GetLength(1) / 2)), 2));
                    if (D > Do)
                    {
                        NewImage.RedImg[u, v] = this.RedImg[u, v];
                        NewImage.GreenImg[u, v] = this.GreenImg[u, v];
                        NewImage.BlueImg[u, v] = this.BlueImg[u, v];

                        NewImage.RedReal[u, v] = this.RedReal[u, v];
                        NewImage.GreenReal[u, v] = this.GreenReal[u, v];
                        NewImage.BlueReal[u, v] = this.BlueReal[u, v];
                    }
                    else
                    {
                        NewImage.RedImg[u, v] = 0;
                        NewImage.GreenImg[u, v] = 0;
                        NewImage.BlueImg[u, v] = 0;

                        NewImage.RedReal[u, v] = 0;
                        NewImage.GreenReal[u, v] = 0;
                        NewImage.BlueReal[u, v] = 0;
                    }
                }
            return NewImage;
        }

        public Image ButterWorthHPF(double Do, double N)
        {
            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);
            NewImage.Fourier = this.Fourier;

            NewImage.RedImg = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];
            NewImage.GreenImg = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];
            NewImage.BlueImg = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];

            NewImage.RedReal = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];
            NewImage.GreenReal = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];
            NewImage.BlueReal = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];

            for (int u = 0; u < this.RedImg.GetLength(0); u++)
                for (int v = 0; v < this.RedImg.GetLength(1); v++)
                {
                    double D = Math.Sqrt(Math.Pow((u - (this.BitmapImage.Width * 1.0 / 2.0)), 2) + Math.Pow((v - (this.BitmapImage.Height * 1.0 / 2.0)), 2));

                    NewImage.RedImg[u, v] = (1.0 - (1.0 / (1.0 + Math.Pow(D / Do, 2.0 * N)))) * this.RedImg[u, v];
                    NewImage.GreenImg[u, v] = (1.0 - (1.0 / (1.0 + Math.Pow(D / Do, 2.0 * N)))) * this.GreenImg[u, v];
                    NewImage.BlueImg[u, v] = (1.0 - (1.0 / (1.0 + Math.Pow(D / Do, 2.0 * N)))) * this.BlueImg[u, v];

                    NewImage.RedReal[u, v] = (1.0 - (1.0 / (1.0 + Math.Pow(D / Do, 2.0 * N)))) * this.RedReal[u, v];
                    NewImage.GreenReal[u, v] = (1.0 - (1.0 / (1.0 + Math.Pow(D / Do, 2.0 * N)))) * this.GreenReal[u, v];
                    NewImage.BlueReal[u, v] = (1.0 - (1.0 / (1.0 + Math.Pow(D / Do, 2.0 * N)))) * this.BlueReal[u, v];
                }

            return NewImage;
        }

        public Image GaussianHPF(double Gamma)
        {
            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);
            NewImage.Fourier = this.Fourier;

            NewImage.RedImg = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];
            NewImage.GreenImg = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];
            NewImage.BlueImg = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];

            NewImage.RedReal = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];
            NewImage.GreenReal = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];
            NewImage.BlueReal = new double[NewImage.Fourier.Width, NewImage.Fourier.Height];


            for (int u = 0; u < this.BitmapImage.Width; u++)
                for (int v = 0; v < this.BitmapImage.Height; v++)
                {
                    double D = Math.Sqrt(Math.Pow((u - (this.BitmapImage.Width * 1.0 / 2.0)), 2) + Math.Pow((v - (this.BitmapImage.Height* 1.0 / 2.0)), 2));

                    double NewR = (1.0 - Math.Exp(-1.0 * Math.Pow(D, 2) / (2.0 * Math.Pow(Gamma, 2)))) * this.RedImg[u, v];
                    double NewG = (1.0 - Math.Exp(-1.0 * Math.Pow(D, 2) / (2.0 * Math.Pow(Gamma, 2)))) * this.GreenImg[u, v];
                    double NewB = (1.0 - Math.Exp(-1.0 * Math.Pow(D, 2) / (2.0 * Math.Pow(Gamma, 2)))) * this.BlueImg[u, v];

                    NewImage.RedImg[u, v] = NewR;
                    NewImage.GreenImg[u, v] = NewG;
                    NewImage.BlueImg[u, v] = NewB;

                    NewR = (1.0 - Math.Exp(-1.0 * Math.Pow(D, 2) / (2.0 * Math.Pow(Gamma, 2)))) * this.RedReal[u, v];
                    NewG = (1.0 - Math.Exp(-1.0 * Math.Pow(D, 2) / (2.0 * Math.Pow(Gamma, 2)))) * this.GreenReal[u, v];
                    NewB = (1.0 - Math.Exp(-1.0 * Math.Pow(D, 2) / (2.0 * Math.Pow(Gamma, 2)))) * this.BlueReal[u, v];

                    NewImage.RedReal[u, v] = NewR;
                    NewImage.GreenReal[u, v] = NewG;
                    NewImage.BlueReal[u, v] = NewB;
                    
                }

            return NewImage;
        }
 
        public Image BandRejectFilter(int Do, int width)
        {
            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);

            int[,] N = new int[BitmapImage.Width, BitmapImage.Height];

            int range1 = Do + (width / 2);
            int range2 = Do - (width / 2);

            for (int u = 0; u < BitmapImage.Width; u++)
            {
                for (int v = 0; v < BitmapImage.Height; v++)
                {
                    double D = Math.Sqrt(Math.Pow((u - (BitmapImage.Width / 2)), 2) + Math.Pow((v - (BitmapImage.Height / 2)), 2));
                    if (D >= range2 && D <= range1)
                    {
                        N[u, v] = 0;
                    }
                    else N[u, v] = 1;
                }
            }

            NewImage.RedReal = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            NewImage.GreenReal = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            NewImage.BlueReal = new double[this.BitmapImage.Width, this.BitmapImage.Height];

            for (int i = 0; i < this.BitmapImage.Width; i++)
            {
                for (int j = 0; j < this.BitmapImage.Height; j++)
                {
                    double NewR, NewG, NewB;

                    NewR = this.RedReal[i, j] * (double)N[i, j];
                    NewG = this.GreenReal[i, j] * (double)N[i, j];
                    NewB = this.BlueReal[i, j] * (double)N[i, j];

                    NewImage.RedReal[i, j] = NewR;
                    NewImage.GreenReal[i, j] = NewG;
                    NewImage.BlueReal[i, j] = NewB;
                }
            }

            NewImage.RedImg = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            NewImage.GreenImg = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            NewImage.BlueImg = new double[this.BitmapImage.Width, this.BitmapImage.Height];

            for (int i = 0; i < this.BitmapImage.Width; i++)
            {
                for (int j = 0; j < this.BitmapImage.Height; j++)
                {
                    double NewR, NewG, NewB;

                    NewR = this.RedImg[i, j] * (double)N[i, j];
                    NewG = this.GreenImg[i, j] * (double)N[i, j];
                    NewB = this.BlueImg[i, j] * (double)N[i, j];

                    NewImage.RedImg[i, j] = NewR;
                    NewImage.GreenImg[i, j] = NewG;
                    NewImage.BlueImg[i, j] = NewB;
                }
            }

            NewImage.Fourier = NewImage.GetMagnitude();
            this.BitmapImage = NewImage.Fourier;
            return NewImage;
        }

        #endregion

        #region Image Restoration

        public Image GeometricMean(int Length)
        {
            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);
            int Border = (Length - 1) / 2;

            for(int i = Border; i < NewImage.BitmapImage.Width - Border; i++)
            {
                for (int j = Border; j < NewImage.BitmapImage.Height - Border; j++)
                {
                    double NewR = 1, NewG = 1, NewB = 1;
                    for (int m = -1 * Border; m <= Border; m++)
                        for (int k = -1 * Border; k <= Border; k++)
                        {
                            NewR *= (BitmapImage.GetPixel(i + m, j + k).R);
                            NewG *= (BitmapImage.GetPixel(i + m, j + k).G);
                            NewB *= (BitmapImage.GetPixel(i + m, j + k).B);
                        }
                    NewR = Math.Pow(NewR, (1.0 / (double)(Length * Length)));
                    NewG = Math.Pow(NewG, (1.0 / (double)(Length * Length)));
                    NewB = Math.Pow(NewB, (1.0 / (double)(Length * Length)));
                    NewImage.BitmapImage.SetPixel(i, j, Color.FromArgb((int)NewR, (int)NewG, (int)NewB));
                }
            }
            IsModified = true;
            return NewImage;
        }

        public Image Median(int Length)
        {
            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);
            BitmapData _bmData = NewImage.BitmapImage.LockBits(new Rectangle(0, 0, this.BitmapImage.Width, this.BitmapImage.Height),ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
            int _stride = _bmData.Stride; 
            IntPtr _Scan0 = _bmData.Scan0;

            bmData = this.BitmapImage.LockBits(new Rectangle(0, 0, this.BitmapImage.Width, this.BitmapImage.Height),
                ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
            stride = bmData.Stride;
            Scan0 = bmData.Scan0;

            int Border = (Length - 1) / 2;
            byte NewR, NewG, NewB;
            unsafe
            {
                for (int i = Border; i < NewImage.BitmapImage.Width - Border; i++)
                {
                    for (int j = Border; j < NewImage.BitmapImage.Height - Border; j++)
                    {
                        ArrayList arrr = new ArrayList();
                        ArrayList arrg = new ArrayList();
                        ArrayList arrb = new ArrayList();

                        for (int m = -1 * Border; m <= Border; m++)
                            for (int k = -1 * Border; k <= Border; k++)
                            {
                                PixelData* pPixel = (PixelData*)(Scan0.ToInt32() + (j + k) * stride + (i + m) * sizeof(PixelData));
                                arrr.Add(pPixel->red);
                                arrg.Add(pPixel->green);
                                arrb.Add(pPixel->blue);
                            }
                        arrr.Sort();
                        arrg.Sort();
                        arrb.Sort();

                        NewR = (byte)arrr[(int)((Length * Length) / 2) + 1];
                        NewG = (byte)arrg[(int)((Length * Length) / 2) + 1];
                        NewB = (byte)arrb[(int)((Length * Length) / 2) + 1];
                        PixelData* _pPixel = (PixelData*)(_Scan0.ToInt32() + (j) * _stride + (i) * sizeof(PixelData));
                        _pPixel->red = NewR;
                        _pPixel->green = NewG;
                        _pPixel->blue = NewB;
                    }
                }
            }
            this.BitmapImage.UnlockBits(bmData);
            NewImage.BitmapImage.UnlockBits(_bmData);
            IsModified = true;
            return NewImage;
        }

        public Image Max(int Length)
        {
            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);
            int Border = (Length - 1) / 2;
            byte NewR, NewG, NewB;

            for (int i = Border; i < NewImage.BitmapImage.Width - Border; i++)
            {
                for (int j = Border; j < NewImage.BitmapImage.Height - Border; j++)
                {
                    ArrayList arrr = new ArrayList();
                    ArrayList arrg = new ArrayList();
                    ArrayList arrb = new ArrayList();

                    for (int m = -1 * Border; m <= Border; m++)
                        for (int k = -1 * Border; k <= Border; k++)
                        {
                            arrr.Add(BitmapImage.GetPixel(i + m, j + k).R);
                            arrg.Add(BitmapImage.GetPixel(i + m, j + k).G);
                            arrb.Add(BitmapImage.GetPixel(i + m, j + k).B);
                        }
                    arrr.Sort();
                    arrg.Sort();
                    arrb.Sort();

                    NewR = (byte)arrr[arrr.Count - 1];
                    NewG = (byte)arrg[arrg.Count - 1];
                    NewB = (byte)arrb[arrb.Count - 1];
                    NewImage.BitmapImage.SetPixel(i, j, Color.FromArgb((int)NewR, (int)NewG, (int)NewB));
                }
            }
            IsModified = true;
            return NewImage;
        }

        public Image Min(int Length)
        {
            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);
            int Border = (Length - 1) / 2;
            byte NewR, NewG, NewB;

            for (int i = Border; i < NewImage.BitmapImage.Width - Border; i++)
            {
                for (int j = Border; j < NewImage.BitmapImage.Height - Border; j++)
                {
                    ArrayList arrr = new ArrayList();
                    ArrayList arrg = new ArrayList();
                    ArrayList arrb = new ArrayList();

                    for (int m = -1 * Border; m <= Border; m++)
                        for (int k = -1 * Border; k <= Border; k++)
                        {
                            arrr.Add(BitmapImage.GetPixel(i + m, j + k).R);
                            arrg.Add(BitmapImage.GetPixel(i + m, j + k).G);
                            arrb.Add(BitmapImage.GetPixel(i + m, j + k).B);
                        }
                    arrr.Sort();
                    arrg.Sort();
                    arrb.Sort();

                    NewR = (byte)arrr[0];
                    NewG = (byte)arrg[0];
                    NewB = (byte)arrb[0];
                    NewImage.BitmapImage.SetPixel(i, j, Color.FromArgb((int)NewR, (int)NewG, (int)NewB));
                }
            }
            IsModified = true;
            return NewImage;
        }

        public Image Midpoint(int Length)
        {
            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);
            int Border = (Length - 1) / 2;
            byte NewR, NewG, NewB;

            for (int i = Border; i < NewImage.BitmapImage.Width - Border; i++)
            {
                for (int j = Border; j < NewImage.BitmapImage.Height - Border; j++)
                {
                    ArrayList arrr = new ArrayList();
                    ArrayList arrg = new ArrayList();
                    ArrayList arrb = new ArrayList();

                    for (int m = -1 * Border; m <= Border; m++)
                        for (int k = -1 * Border; k <= Border; k++)
                        {
                            arrr.Add(BitmapImage.GetPixel(i + m, j + k).R);
                            arrg.Add(BitmapImage.GetPixel(i + m, j + k).G);
                            arrb.Add(BitmapImage.GetPixel(i + m, j + k).B);
                        }
                    arrr.Sort();
                    arrg.Sort();
                    arrb.Sort();

                    NewR = (byte)(((byte)arrr[0] + (byte)arrr[arrr.Count - 1]) /2);
                    NewG = (byte)(((byte)arrg[0] + (byte)arrg[arrg.Count - 1]) / 2);
                    NewB = (byte)(((byte)arrb[0] + (byte)arrb[arrb.Count - 1]) /2);
                    NewImage.BitmapImage.SetPixel(i, j, Color.FromArgb((int)NewR, (int)NewG, (int)NewB));
                }
            }
            IsModified = true;
            return NewImage;
        }

        public Image AdaptiveMedian(int MaxSize)
        {
            int CenterPixR, CenterPixG, CenterPixB;
            int size = 3;
            int[,] Mask = new int[size, size];
            byte MaxR, MaxG, MaxB, MinR, MinG, MinB, MedR, MedG, MedB;

            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);

            byte NewR = 0, NewG = 0, NewB = 0;

            for (int i = 0; i < BitmapImage.Width; i++)
            {
                for (int j = 0; j < BitmapImage.Height; j++)
                {
                    size = 3;

                Sorting:
                    ArrayList arrr = new ArrayList();
                    ArrayList arrg = new ArrayList();
                    ArrayList arrb = new ArrayList();

                    int Border = (size) / 2;
                    for (int m = -1 * Border; m <= Border; m++)
                    {
                        for (int k = -1 * Border; k <= Border; k++)
                        {
                            if ((m + i) >= 0 && (m + i) < BitmapImage.Width && (k + j) >= 0 && (k + j) < BitmapImage.Height)
                            {
                                Color c = BitmapImage.GetPixel(m + i, k + j);
                                arrr.Add(c.R);
                                arrg.Add(c.G);
                                arrb.Add(c.B);
                            }
                        }
                    }
                    CenterPixR = (byte)arrr[(arrr.Count - 1) / 2];
                    CenterPixG = (byte)arrg[(arrg.Count - 1) / 2];
                    CenterPixB = (byte)arrb[(arrb.Count - 1) / 2];

                    arrr.Sort();
                    arrg.Sort();
                    arrb.Sort();

                    MinR = (byte)arrr[0];
                    MaxR = (byte)arrr[arrr.Count - 1];
                    MedR = (byte)arrr[(int)((arrr.Count - 1) / 2)];

                    MinG = (byte)arrg[0];
                    MaxG = (byte)arrg[arrg.Count - 1];
                    MedG = (byte)arrg[(int)((arrg.Count - 1) / 2)];

                    MinB = (byte)arrb[0];
                    MaxB = (byte)arrb[arrb.Count - 1];
                    MedB = (byte)arrb[(int)((arrb.Count - 1) / 2)];

                    if (MinR < MedR && MedR < MaxR)
                    {
                        if (MinG < MedG && MedG < MaxG)
                        {
                            if (MinB < MedB && MedB < MaxB)
                            {
                                //Color c = BitmapImage.GetPixel(i, j);
                                if (MinR < CenterPixR && CenterPixR < MaxR)
                                {
                                    NewR = (byte)CenterPixR;//c.R;//CenterPix = CenterPix;
                                }
                                else
                                {
                                    NewR = (byte)MedR;
                                }
                                if (MinG < CenterPixG && CenterPixG < MaxG)
                                {
                                    NewG = (byte)CenterPixG;//c.G;//CenterPix = CenterPix;                        
                                }
                                else
                                {
                                    NewG = (byte)MedG;
                                }
                                if (MinB < CenterPixB && CenterPixB < MaxB)
                                {
                                    NewB = (byte)CenterPixB;//c.B;//CenterPix = CenterPix;
                                }
                                else
                                {
                                    NewB = (byte)MedB;
                                }
                            }
                        }
                    }
                    else
                    {
                        size += 2;
                        if (size <= MaxSize)
                        {
                            goto Sorting;
                        }
                        else
                        {
                            size -= 2;
                            //Color c = BitmapImage.GetPixel(i, j);
                            NewR = (byte)CenterPixR;//c.R;
                            NewG = (byte)CenterPixG;//c.G;
                            NewB = (byte)CenterPixB;//c.B;
                        }//else CenterPix = CenterPix;
                    }

                    NewImage.BitmapImage.SetPixel(i, j, Color.FromArgb((int)NewR, (int)NewG, (int)NewB));
                }
            }

            return NewImage;
        }

        public Image Alphatrim(int Alpha, int MaskSize)
        {
            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(BitmapImage.Width, BitmapImage.Height);

            int Border = (MaskSize - 1) / 2;
            //byte NewR, NewG, NewB;

            if ((MaskSize % 2) == 0)
            {
                MaskSize += 1;
            }

            //double[,] BufferR = new double[BitmapImage.Width, BitmapImage.Height];
            //double[,] BufferB = new double[BitmapImage.Width, BitmapImage.Height];
            //double[,] BufferG = new double[BitmapImage.Width, BitmapImage.Height];

            for (int i = Border; i < NewImage.BitmapImage.Width - Border; i++)//if masksize=3 and border=1
            {
                for (int j = Border; j < NewImage.BitmapImage.Height - Border; j++)
                {
                    ArrayList arrr = new ArrayList();
                    ArrayList arrg = new ArrayList();
                    ArrayList arrb = new ArrayList();

                    byte Avgr = 0, Avgg = 0, Avgb = 0;

                    for (int m = -Border; m <= Border; m++)
                    {
                        for (int k = -Border; k <= Border; k++)
                        {
                            arrr.Add(BitmapImage.GetPixel(i + m, j + k).R);
                            arrg.Add(BitmapImage.GetPixel(i + m, j + k).G);
                            arrb.Add(BitmapImage.GetPixel(i + m, j + k).B);
                        }
                    }
                    arrr.Sort();
                    arrg.Sort();
                    arrb.Sort();

                    arrr.RemoveRange(0, Alpha - 1); arrr.RemoveRange(arrr.Count - Alpha - 1, Alpha);
                    arrg.RemoveRange(0, Alpha - 1); arrg.RemoveRange(arrr.Count - Alpha - 1, Alpha);
                    arrb.RemoveRange(0, Alpha - 1); arrb.RemoveRange(arrr.Count - Alpha - 1, Alpha);

                    for (int a = 0; a < arrr.Count; a++)//Red
                        Avgr += (byte)arrr[a];
                    Avgr /= (byte)arrr.Count;

                    for (int a = 0; a < arrg.Count; a++)//Green
                        Avgg += (byte)arrg[a];
                    Avgg /= (byte)arrg.Count;

                    for (int a = 0; a < arrb.Count; a++)//Blue
                        Avgb += (byte)arrb[a];
                    Avgb /= (byte)arrb.Count;

                    NewImage.BitmapImage.SetPixel(i, j, Color.FromArgb((int)Avgr, (int)Avgg, (int)Avgb));
                    //BufferR[i, j] = Avgr;
                    //BufferG[i, j] = Avgg;
                    //BufferB[i, j] = Avgb;
                }
            }
            //NewImage.BitmapImage = Normalize(BufferR, BufferB, BufferG, 0, 255);
            return NewImage;
        }

        #endregion

        #region Noise

        public Image AdditiveNoise(string NoiseType, int v1, int v2, int percentage)
        {
            double NumOfPix;
            Random RandNum = new Random();
            int W = BitmapImage.Width;
            int H = BitmapImage.Height;
            double[,] Noise = new double[W, H];

            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap((int)BitmapImage.Width, (int)BitmapImage.Height);

            if (NoiseType == "Uniform")
            {
                NumOfPix = (1.0 / (v2 - v1)) * BitmapImage.Width * BitmapImage.Height * (percentage / 100.0);

                for (int i = v1; i <= v2; i++) // from a to b
                {
                    for (int j = 0; j < (int)NumOfPix; j++)
                    {
                        int x = RandNum.Next(0, BitmapImage.Width - 1);
                        int y = RandNum.Next(0, BitmapImage.Height - 1);
                        Noise[x, y] = i;
                    }
                }
            }
            else if (NoiseType == "Gaussian")
            {
                for (int i = 0; i < 255; i++) // v1 = m , v2 = Stand Dev.
                {
                    NumOfPix = (1.0 / ((Math.Sqrt(2 * Math.PI)) * v2)) * (Math.Exp((-Math.Pow((i - v1), 2)) / (2 * Math.Pow(v2, 2)))) * BitmapImage.Width * BitmapImage.Height * (percentage / 100.0);
                    for (int j = 0; j < (int)NumOfPix; j++)
                    {
                        int x = RandNum.Next(0, BitmapImage.Width);
                        int y = RandNum.Next(0, BitmapImage.Height);
                        Noise[x, y] = i;
                    }
                }
            }

            PreparingNormalization(Noise);

            for (int X = 0; X < NewImage.BitmapImage.Width; X++)
            {
                for (int Y = 0; Y < NewImage.BitmapImage.Height; Y++)
                {
                    Color OldPixel = this.BitmapImage.GetPixel(X, Y);
                    int NoisePixel = (int)Noise[X, Y];

                    int R = OldPixel.R + NoisePixel;
                    int G = OldPixel.G + NoisePixel;
                    int B = OldPixel.B + NoisePixel;

                    //Normalization
                    R = Normalize(R, MinR, MaxR, 0, 255);
                    G = Normalize(G, MinG, MaxG, 0, 255);
                    B = Normalize(B, MinB, MaxB, 0, 255);

                    try { NewImage.BitmapImage.SetPixel(X, Y, Color.FromArgb(R, G, B)); }
                    catch { }
                }
            }

            return NewImage;
        }

        public Image PeriodicNoise(int A, int Uo, int Vo, int Bx, int By)
        {
            int W = BitmapImage.Width;
            int H = BitmapImage.Height;
            double[,] Noise = new double[W, H];

            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(BitmapImage.Width, BitmapImage.Height);

            for (int x = 0; x < W; x++)
            {
                for (int y = 0; y < H; y++)
                {
                    double firstpart = ((2 * Math.PI * Uo) * (x + Bx)) / W;
                    double secondpart = ((2 * Math.PI * Vo) * (y + By)) / H;
                    double NoiseVal = A * (Math.Sin(firstpart + secondpart));
                    Noise[x, y] = NoiseVal;
                }
            }

            PreparingNormalization(Noise);

            for (int X = 0; X < NewImage.BitmapImage.Width; X++)
            {
                for (int Y = 0; Y < NewImage.BitmapImage.Height; Y++)
                {
                    Color OldPixel = this.BitmapImage.GetPixel(X, Y);
                    double NoisePixel = Noise[X, Y];

                    double R = OldPixel.R + NoisePixel;
                    double G = OldPixel.G + NoisePixel;
                    double B = OldPixel.B + NoisePixel;

                    //Normalization
                    R = Normalize((int)R, MinR, MaxR, 0, 255);
                    G = Normalize((int)G, MinG, MaxG, 0, 255);
                    B = Normalize((int)B, MinB, MaxB, 0, 255);

                    try { NewImage.BitmapImage.SetPixel(X, Y, Color.FromArgb((int)R, (int)G, (int)B)); }
                    catch { }
                }
            }

            return NewImage;
        }
        public Image SaltPepper(float Ps, float Pp)
        {
            Image NewImage = new Image();
            NewImage.BitmapImage = this.BitmapImage;

            float NsFloat = Ps * this.BitmapImage.Width * this.BitmapImage.Height;
            float NpFloat = Pp * this.BitmapImage.Width * this.BitmapImage.Height;

            int Ns = (int)NsFloat; int Np = (int)NpFloat;

            Random random = new Random(1);

            int RandomPixelX; int RandomPixelY;

            for (int i = 0; i < Ns; i++)
            {
                RandomPixelX = random.Next(0, NewImage.BitmapImage.Width - 1);
                RandomPixelY = random.Next(0, NewImage.BitmapImage.Height - 1);

                NewImage.BitmapImage.SetPixel(RandomPixelX, RandomPixelY, Color.FromArgb(255, 255, 255));
            }
            for (int i = 0; i < Np; i++)
            {
                RandomPixelX = random.Next(0, NewImage.BitmapImage.Width - 1);
                RandomPixelY = random.Next(0, NewImage.BitmapImage.Height - 1);

                NewImage.BitmapImage.SetPixel(RandomPixelX, RandomPixelY, Color.FromArgb(0, 0, 0));
            }
            return NewImage;
        }

        #endregion Noise

        #region Using MATLAB

        public Image Multiplying2ImagesUsingMatlab(Image img1, Image img2)
        {
            MatlabMultiplication.MatlabMultiplicationclassClass m = new MatlabMultiplication.MatlabMultiplicationclassClass();
            img2 = img2.ImageResize(img1.BitmapImage.Width, img1.BitmapImage.Height);
            byte[, ,] I1 = new byte[img1.BitmapImage.Height + 1, img1.BitmapImage.Width + 1, 3];
            byte[, ,] I2 = new byte[img1.BitmapImage.Height + 1, img1.BitmapImage.Width + 1, 3];
            byte[, ,] I3 = new byte[img1.BitmapImage.Height + 1, img1.BitmapImage.Width + 1, 3];

            for (int j = 1; j <= img1.BitmapImage.Height; j++)
                for (int i = 1; i <= img1.BitmapImage.Width; i++)
                {
                    I1[j, i, 0] = img1.BitmapImage.GetPixel(i - 1, j - 1).R;
                    I1[j, i, 1] = img1.BitmapImage.GetPixel(i - 1, j - 1).G;
                    I1[j, i, 2] = img1.BitmapImage.GetPixel(i - 1, j - 1).B;

                    I2[j, i, 0] = img2.BitmapImage.GetPixel(i - 1, j - 1).R;
                    I2[j, i, 1] = img2.BitmapImage.GetPixel(i - 1, j - 1).G;
                    I2[j, i, 2] = img2.BitmapImage.GetPixel(i - 1, j - 1).B;
                }

            Object obj = (Object)I3;
            m.multiplytwoimages(1, ref obj, I1, I2);

            I3 = (byte[, ,])obj;
            Bitmap bit = new Bitmap(img1.BitmapImage.Width, img2.BitmapImage.Height);

            for (int j = 1; j <= img1.BitmapImage.Height; j++)
                for (int i = 1; i <= img1.BitmapImage.Width; i++)
                {
                    bit.SetPixel(i - 1, j - 1, Color.FromArgb((int)I3[j, i, 1], (int)I3[j, i, 2], (int)I3[j, i, 3]));
                }

            Image NewImage = new Image(bit);
            return NewImage;
        }

        public Image ApplyFourier()
        {
            double[,] Rr = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            double[,] Gg = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            double[,] Bb = new double[this.BitmapImage.Width, this.BitmapImage.Height];

            RedReal = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            GreenReal = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            BlueReal = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            RedImg = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            GreenImg = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            BlueImg = new double[this.BitmapImage.Width, this.BitmapImage.Height];

            for (int i = 0; i < this.BitmapImage.Height; i++)
            {
                for (int j = 0; j < this.BitmapImage.Width; j++)
                {
                    Rr[j, i] = (double)this.BitmapImage.GetPixel(j, i).R;
                    Gg[j, i] = (double)this.BitmapImage.GetPixel(j, i).G;
                    Bb[j, i] = (double)this.BitmapImage.GetPixel(j, i).B;
                }
            }

            EngMATAccess engmat = new EngMATAccess();
            engmat.SetMatrix("R", Rr);
            engmat.SetMatrix("G", Gg);
            engmat.SetMatrix("B", Bb); engmat.Evaluate("RedReal=real(fftshift(fft2(R)));GreenReal=real(fftshift(fft2(G)));BlueReal=real(fftshift(fft2(B)));");
            engmat.GetMatrix("RedReal", ref RedReal);
            engmat.GetMatrix("GreenReal", ref GreenReal);
            engmat.GetMatrix("BlueReal", ref BlueReal);

            engmat.Evaluate("RedImg=imag(fftshift(fft2(R)));GreenImg=imag(fftshift(fft2(G)));BlueImg=imag(fftshift(fft2(B)));");
            engmat.GetMatrix("RedImg", ref RedImg);
            engmat.GetMatrix("GreenImg", ref GreenImg);
            engmat.GetMatrix("BlueImg", ref BlueImg);

            Fourier = GetMagnitude() ;
            return this;
        }

        public Bitmap GetMagnitude()
        {
            MagR = new double[BitmapImage.Width, BitmapImage.Height];
            MagG = new double[BitmapImage.Width, BitmapImage.Height];        
            MagB = new double[BitmapImage.Width, BitmapImage.Height];

            Bitmap MagBm = new Bitmap(BitmapImage.Width, BitmapImage.Height);

            for (int i = 0; i < BitmapImage.Height; i++)
            {
                for (int j = 0; j < BitmapImage.Width; j++)
                {
                    MagR[j, i] = Math.Sqrt((RedReal[j, i] * RedReal[j, i]) + (RedImg[j, i] * RedImg[j, i]));
                    MagG[j, i] = Math.Sqrt((GreenReal[j, i] * GreenReal[j, i]) + (GreenImg[j, i] * GreenImg[j, i]));
                    MagB[j, i] = Math.Sqrt((BlueReal[j, i] * BlueReal[j, i]) + (BlueImg[j, i] * BlueImg[j, i]));
                }
            }

            for (int i = 0; i < BitmapImage.Height; i++)
            {
                for (int j = 0; j < BitmapImage.Width; j++)
                {
                    MagR[j, i] = Math.Log(MagR[j, i] + 1);
                    MagG[j, i] = Math.Log(MagG[j, i] + 1);
                    MagB[j, i] = Math.Log(MagB[j, i] + 1);
                }
            }
            return Normalize(MagR,MagG,MagB,0,255);
        }

        public Image ApplyFourierInverse()
        {
            double[,] InvRed = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            double[,] InvGreen = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            double[,] InvBlue = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            Color[,] InvBuffer = new Color[this.BitmapImage.Width, this.BitmapImage.Height];

            EngMATAccess engmat = new EngMATAccess();
            engmat.SetMatrix("RedReal", RedReal);
            engmat.SetMatrix("GreenReal", GreenReal);
            engmat.SetMatrix("BlueReal", BlueReal);
            engmat.SetMatrix("RedImg", RedImg);
            engmat.SetMatrix("GreenImg", GreenImg);
            engmat.SetMatrix("BlueImg", BlueImg);
            engmat.Evaluate("InvRed=real(ifft2(ifftshift(RedReal+i*RedImg)));InvGreen=real(ifft2(ifftshift(GreenReal+i*GreenImg)));InvBlue=real(ifft2(ifftshift(BlueReal+i*BlueImg)));");
            engmat.GetMatrix("InvRed", ref InvRed);
            engmat.GetMatrix("InvGreen", ref InvGreen);
            engmat.GetMatrix("InvBlue", ref InvBlue);

            BitmapImage = Normalize(InvRed, InvGreen, InvBlue, 0, 255);

            engmat.Close();
            return this;
        }

        #endregion Using MATLAB

        #region Helper Functions

        public Color CutOff(int R, int G, int B)
        {
            int Min = 0; int Max = 255;

            if (R < Min) R = Min;
            else if (R > Max) R = Max;
            if (G < Min) G = Min;
            else if (G > Max) G = Max;
            if (B < Min) B = Min;
            else if (B > Max) B = Max;

            return Color.FromArgb(R, G, B);
        }

        public void PreparingNormalization(double[,] Noisy)
        {            
            MinR = BitmapImage.GetPixel(0, 0).R + (int)Noisy[0, 0];
            MaxR = BitmapImage.GetPixel(0, 0).R + (int)Noisy[0, 0];
            MinG = BitmapImage.GetPixel(0, 0).G + (int)Noisy[0, 0];
            MaxG = BitmapImage.GetPixel(0, 0).G + (int)Noisy[0, 0];
            MinB = BitmapImage.GetPixel(0, 0).B + (int)Noisy[0, 0];
            MaxB = BitmapImage.GetPixel(0, 0).B + (int)Noisy[0, 0];

            // Getting The Min & Max Value in The New Bitmap
            for (int i = 0; i < BitmapImage.Width; i++)
            {
                for (int j = 0; j < BitmapImage.Height; j++)
                {
                    Color OldPixel = this.BitmapImage.GetPixel(i, j);
                    int NoisePixel = (int)Noisy[i, j];

                    int R = OldPixel.R + NoisePixel;
                    int G = OldPixel.G + NoisePixel;
                    int B = OldPixel.B + NoisePixel;

                    if (R > MaxR) MaxR = R;
                    else if (R < MinR) MinR = R;
                    if (G > MaxG) MaxG = G;
                    else if (G < MinG) MinG = G;
                    if (B > MaxB) MaxB = B;
                    else if (B < MinB) MinB = B;
                }
            }
        }

        public int Normalize(int OldVal, int OldMin, int OldMax, int NewMin, int NewMax)
        {
            return (int)((float)((float)(OldVal - OldMin) / (float)(OldMax - OldMin)) * (NewMax - NewMin) + NewMin);
        }

        public Bitmap Normalize(double[,] OldValR, double[,] OldValG, double[,] OldValB, int NewMin, int NewMax)
        {
            Bitmap New = new Bitmap(OldValR.GetLength(0), OldValR.GetLength(1));
            double OldMinR = GetMinVal(OldValR); double OldMaxR = GetMaxVal(OldValR);
            double OldMinG = GetMinVal(OldValG); double OldMaxG = GetMaxVal(OldValG);
            double OldMinB = GetMinVal(OldValB); double OldMaxB = GetMaxVal(OldValB);

            int NewR; int NewG; int NewB;
            
            for (int i = 0; i < New.Width; i++)
                for (int j = 0; j < New.Height; j++)
                {
                    NewR = (int)((float)((float)(OldValR[i, j] - OldMinR) / (float)(OldMaxR - OldMinR)) * (NewMax - NewMin) + NewMin);
                    NewG = (int)((float)((float)(OldValG[i, j] - OldMinG) / (float)(OldMaxG - OldMinG)) * (NewMax - NewMin) + NewMin);
                    NewB = (int)((float)((float)(OldValB[i, j] - OldMinB) / (float)(OldMaxB - OldMinB)) * (NewMax - NewMin) + NewMin);

                    //In case of power, after normalization we need cutoff
                    if (NewR > 255 || NewR < 0) NewR = 255;
                    if (NewG > 255 || NewG < 0) NewG = 255;
                    if (NewB > 255 || NewB < 0) NewB = 255;

                    New.SetPixel(i, j, Color.FromArgb(NewR, NewG, NewB));
                }

            return New;
        }

        public double[,] Normalize(double[,] OldVal, int NewMin, int NewMax)
        {
            double[,] New = new double[OldVal.GetLength(0), OldVal.GetLength(1)];
            double OldMin = GetMinVal(OldVal); double OldMax = GetMaxVal(OldVal);
            int NewPix;
            for (int i = 0; i < New.GetLength(0); i++)
                for (int j = 0; j < New.GetLength(1); j++)
                {
                    NewPix = (int)((float)((float)(OldVal[i, j] - OldMin) / (float)(OldMax - OldMin)) * (NewMax - NewMin) + NewMin);

                    //In case of power, after normalization we need cutoff
                    if (NewPix > 255 || NewPix < 0) NewPix = 255;

                    New[i, j] = NewPix;
                }

            return New;

        }

        public void Normalizedouble(double[,] R, double[,] G, double[,] B)
        {
            double Old_minR = 99999999, Old_maxR = -1;
            double Old_minG = 99999999, Old_maxG = -1;
            double Old_minB = 99999999, Old_maxB = -1;

            for (int i = 0; i < BitmapImage.Width; i++)
            {
                for (int j = 0; j < BitmapImage.Height; j++)
                {
                    if (R[i, j] < Old_minR)
                        Old_minR = R[i, j];
                    if (R[i, j] > Old_maxR)
                        Old_maxR = R[i, j];

                    if (G[i, j] < Old_minG)
                        Old_minG = G[i, j];
                    if (G[i, j] > Old_maxG)
                        Old_maxG = G[i, j];

                    if (B[i, j] < Old_minB)
                        Old_minB = B[i, j];
                    if (B[i, j] > Old_maxB)
                        Old_maxB = B[i, j];

                }
            }

            for (int i = 0; i < BitmapImage.Width; i++)
            {
                for (int j = 0; j < BitmapImage.Height; j++)
                {

                    R[i, j] = ((R[i, j] - Old_minR) / ((Old_maxR - Old_minR) + 1)) * (255 - 0) + 0;
                    G[i, j] = ((G[i, j] - Old_minG) / ((Old_maxG - Old_minG) + 1)) * (255 - 0) + 0;
                    B[i, j] = ((B[i, j] - Old_minB) / ((Old_maxB - Old_minB) + 1)) * (255 - 0) + 0;

                }
            }

        }

        public int AllOperations(int R1, int R2, string Operation)
        {
            int R;

            switch (Operation)
            {
                case "Addition":
                    R = R1 + R2;
                    return R;

                case "Subtraction":
                    R = R1 - R2;
                    return R;

                case "Multiplication":
                    R = R1 * R2;
                    return R;

                case "Division":
                    R = R1 / R2;
                    return R;
                case "Or":
                    R = R1 | R2;
                    return R;
                case "And":
                    R = R1 & R2;
                    return R;
                default:
                    return 0;
            }
        }

        public Color GetMaxVal(Bitmap Map)
        {
            int MaxR = Map.GetPixel(0, 0).R;
            int MaxG = Map.GetPixel(0, 0).G;
            int MaxB = Map.GetPixel(0, 0).B; 

            for (int i = 0; i < Map.Width; i++)
                for (int j = 0; j < Map.Height; j++)
                {
                    Color Pixel = Map.GetPixel(i, j);

                    if (Pixel.R > MaxR) MaxR = Pixel.R;
                    if (Pixel.G > MaxG) MaxG = Pixel.G;
                    if (Pixel.B > MaxB) MaxB = Pixel.B;
                }
            return Color.FromArgb(MaxR, MaxG, MaxB);
        }

        public int GetMinVal(int First, int Sec)
        {
            if (First < Sec) return First;
            else return Sec;
        }

        public Color GetMinVal(Bitmap Map)
        {
            int MinR = Map.GetPixel(0, 0).R;
            int MinG = Map.GetPixel(0, 0).G;
            int MinB = Map.GetPixel(0, 0).B;

            for (int i = 0; i < Map.Width; i++)
                for (int j = 0; j < Map.Height; j++)
                {
                    Color Pixel = Map.GetPixel(i, j);

                    if (Pixel.R < MinR) MinR = Pixel.R;
                    if (Pixel.G < MinG) MinG = Pixel.G;
                    if (Pixel.B < MinB) MinB = Pixel.B;
                }
            return Color.FromArgb(MinR, MinG, MinB);
        }

        public double GetMaxVal(double[,] Arr)
        {
            double Max = 0;
            for (int i = 0; i < Arr.GetLength(0); i++)
                for (int j = 0; j < Arr.GetLength(1); j++)
                    if (Arr[i, j] > Max) Max = Arr[i, j];

            return Max;
        }

        public double GetMinVal(double[,] Arr)
        {
            double Min = double.MaxValue;
            for (int i = 0; i < Arr.GetLength(0); i++)
                for (int j = 0; j < Arr.GetLength(1); j++)
                {
                    if (Arr[i, j] < Min)
                    {
                        Min = Arr[i, j];
                    }
                }
            return Min;
        }

        #endregion Helper Functions

        #region Extra

        public Bitmap Crop(Point p1, Point p2)
        {
            int NewWidth = p2.X - p1.X;
            int NewHeight = p2.Y - p1.Y;
            int starti = p1.X,startj = p1.Y;
            if (NewWidth > 0 && NewHeight < 0)
            {
                starti = p2.X;
                startj = p2.Y;
            }
            else if (NewWidth < 0 && NewHeight > 0)
            {
                starti = p2.X;
                startj = p1.Y;
            }
            else if (NewWidth < 0 && NewHeight < 0)
            {
                starti = p1.X;
                startj = p2.Y;
            }

            Bitmap NewImage = new Bitmap(NewWidth, NewHeight);

            for (int i = starti; i < starti + NewWidth; i++)
            {
                for (int j = startj; j < startj + NewHeight; j++)
                {
                    NewImage.SetPixel(i - starti, j - startj, this.BitmapImage.GetPixel(i, j));
                }
            }
            return NewImage;
        }

        public Bitmap GrayScale(Bitmap Old)
        {
            Bitmap New = new Bitmap(Old.Width, Old.Height);

            for (int i = 0; i < New.Width; i++)
                for (int j = 0; j < New.Height; j++)
                {
                    int avg = (Old.GetPixel(i, j).R + Old.GetPixel(i, j).G + Old.GetPixel(i, j).B) / 3;
                    New.SetPixel(i, j, Color.FromArgb(avg, avg, avg));
                }

            return New;
        }

        public Image ImageLogic(Image First, Image Second, string Operation)
        {
            // Resizing Second  to First
            Second.ImageResize(First.BitmapImage.Width, First.BitmapImage.Height);
            Image NewImage = new Image(new Bitmap(First.BitmapImage.Width, First.BitmapImage.Height));
            int OldMinR = AllOperations(First.BitmapImage.GetPixel(0, 0).R, Second.BitmapImage.GetPixel(0, 0).R, Operation);
            int OldMaxR = AllOperations(First.BitmapImage.GetPixel(0, 0).R, Second.BitmapImage.GetPixel(0, 0).R, Operation);
            int OldMinG = AllOperations(First.BitmapImage.GetPixel(0, 0).G, Second.BitmapImage.GetPixel(0, 0).G, Operation);
            int OldMaxG = AllOperations(First.BitmapImage.GetPixel(0, 0).G, Second.BitmapImage.GetPixel(0, 0).G, Operation);
            int OldMinB = AllOperations(First.BitmapImage.GetPixel(0, 0).B, Second.BitmapImage.GetPixel(0, 0).B, Operation);
            int OldMaxB = AllOperations(First.BitmapImage.GetPixel(0, 0).B, Second.BitmapImage.GetPixel(0, 0).B, Operation);

            // Getting The Min & Max Value in The New Bitmap
            for (int i = 0; i < NewImage.BitmapImage.Width; i++)
                for (int j = 0; j < NewImage.BitmapImage.Height; j++)
                {
                    Color OldPixel1 = First.BitmapImage.GetPixel(i, j);
                    Color OldPixel2 = Second.BitmapImage.GetPixel(i, j);

                    //int A = AllOperations(OldPixel1.A, OldPixel2.A, Operation);
                    int R = AllOperations(OldPixel1.R, OldPixel2.R, Operation);
                    int G = AllOperations(OldPixel1.G, OldPixel2.G, Operation);
                    int B = AllOperations(OldPixel1.B, OldPixel2.B, Operation);

                    //if (A > OldMaxA) OldMaxA = A;
                    //else if (A < OldMinA) OldMinA = A;
                    if (R > OldMaxR) OldMaxR = R;
                    else if (R < OldMinR) OldMinR = R;
                    if (G > OldMaxG) OldMaxG = G;
                    else if (G < OldMinG) OldMinG = G;
                    if (B > OldMaxB) OldMaxB = B;
                    else if (B < OldMinB) OldMinB = B;
                }

            // Normalizing The New Bitmap
            for (int i = 0; i < NewImage.BitmapImage.Width; i++)
                for (int j = 0; j < NewImage.BitmapImage.Height; j++)
                {
                    Color OldPixel1 = First.BitmapImage.GetPixel(i, j);
                    Color OldPixel2 = Second.BitmapImage.GetPixel(i, j);

                    int R = System.Math.Abs(AllOperations(OldPixel1.R, OldPixel2.R, Operation));
                    int G = System.Math.Abs(AllOperations(OldPixel1.G, OldPixel2.G, Operation));
                    int B = System.Math.Abs(AllOperations(OldPixel1.B, OldPixel2.B, Operation));

                    R = Normalize(R, OldMinR, OldMaxR, 0, 255);
                    G = Normalize(G, OldMinG, OldMaxG, 0, 255);
                    B = Normalize(B, OldMinB, OldMaxB, 0, 255);

                    try { NewImage.BitmapImage.SetPixel(i, j, Color.FromArgb(R, G, B)); }
                    catch { }
                }
            this.BitmapImage = NewImage.BitmapImage;
            return NewImage;
        }

        public Image BandPassFilter(int Do, int width)
        {
            Image NewImage = new Image();
            NewImage.BitmapImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);

            int[,] N = new int[BitmapImage.Width, BitmapImage.Height];

            int range1 = Do + (width / 2);
            int range2 = Do - (width / 2);

            for (int u = 0; u < BitmapImage.Width; u++)
            {
                for (int v = 0; v < BitmapImage.Height; v++)
                {
                    double D = Math.Sqrt(Math.Pow((u - (BitmapImage.Width / 2)), 2) + Math.Pow((v - (BitmapImage.Height / 2)), 2));
                    if (D >= range2 && D <= range1)
                    {
                        N[u, v] = 1;
                    }
                    else N[u, v] = 0;
                }
            }

            NewImage.RedReal = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            NewImage.GreenReal = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            NewImage.BlueReal = new double[this.BitmapImage.Width, this.BitmapImage.Height];

            for (int i = 0; i < this.BitmapImage.Width; i++)
            {
                for (int j = 0; j < this.BitmapImage.Height; j++)
                {
                    double NewR, NewG, NewB;

                    NewR = this.RedReal[i, j] * (double)N[i, j];
                    NewG = this.GreenReal[i, j] * (double)N[i, j];
                    NewB = this.BlueReal[i, j] * (double)N[i, j];

                    NewImage.RedReal[i, j] = NewR;
                    NewImage.GreenReal[i, j] = NewG;
                    NewImage.BlueReal[i, j] = NewB;
                }
            }

            NewImage.RedImg = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            NewImage.GreenImg = new double[this.BitmapImage.Width, this.BitmapImage.Height];
            NewImage.BlueImg = new double[this.BitmapImage.Width, this.BitmapImage.Height];

            for (int i = 0; i < this.BitmapImage.Width; i++)
            {
                for (int j = 0; j < this.BitmapImage.Height; j++)
                {
                    double NewR, NewG, NewB;

                    NewR = this.RedImg[i, j] * (double)N[i, j];
                    NewG = this.GreenImg[i, j] * (double)N[i, j];
                    NewB = this.BlueImg[i, j] * (double)N[i, j];

                    NewImage.RedImg[i, j] = NewR;
                    NewImage.GreenImg[i, j] = NewG;
                    NewImage.BlueImg[i, j] = NewB;
                }
            }

            NewImage.Fourier = NewImage.GetMagnitude();
            this.BitmapImage = NewImage.Fourier;
            return NewImage;
        }

        public void NotchReject(int x, int y, int NotchWidth)
        {
            int x2 = this.BitmapImage.Width - x - 2 , y2 = this.BitmapImage.Height - y - 2;

            for (int i = 0; i < this.BitmapImage.Width; i++)
            {
                for (int j = 0; j < this.BitmapImage.Height; j++)
                {
                    double D1 = Math.Sqrt(Math.Pow(i - x,2) + Math.Pow(j - y,2));
                    double D2 = Math.Sqrt(Math.Pow(i - x2,2) + Math.Pow(j - y2,2));

                    if ((D1 < (double)(NotchWidth*1.0 / 1.0) || D2 < (double)(NotchWidth*1.0 / 1.0)))
                    {
                        RedReal[i, j] = 0;
                        GreenReal[i, j] = 0;
                        BlueReal[i, j] = 0;

                        RedImg[i, j] = 0;
                        GreenImg[i, j] = 0;
                        BlueImg[i, j] = 0;
                    }
                }
            }
        }

        public void NotchPass(int x, int y, int NotchWidth)
        {
            int x2 = this.BitmapImage.Width - x - 2, y2 = this.BitmapImage.Height - y - 2;

            for (int i = 0; i < this.BitmapImage.Width; i++)
            {
                for (int j = 0; j < this.BitmapImage.Height; j++)
                {
                    double D1 = Math.Sqrt(Math.Pow(i - x, 2) + Math.Pow(j - y, 2));
                    double D2 = Math.Sqrt(Math.Pow(i - x2, 2) + Math.Pow(j - y2, 2));

                    if ((D1 < (double)(NotchWidth * 1.0 / 1.0) || D2 < (double)(NotchWidth * 1.0 / 1.0)))
                    {
                    }
                    else
                    {
                        RedReal[i, j] = 0;
                        GreenReal[i, j] = 0;
                        BlueReal[i, j] = 0;

                        RedImg[i, j] = 0;
                        GreenImg[i, j] = 0;
                        BlueImg[i, j] = 0;
                    }
                }
            }
        }

        public Image BlackWhite()
        {
            Bitmap NewImage = new Bitmap(this.BitmapImage.Width, this.BitmapImage.Height);

            for (int i = 0; i < this.BitmapImage.Width; i++)
            {
                for (int j = 0; j < this.BitmapImage.Height; j++)
                {
                    if ((this.BitmapImage.GetPixel(i, j).R + this.BitmapImage.GetPixel(i, j).G + this.BitmapImage.GetPixel(i, j).B)/3.0 < 125)
                        NewImage.SetPixel(i, j, Color.FromArgb(0,0,0));
                    else
                        NewImage.SetPixel(i, j, Color.FromArgb(255, 255, 255));
                }
            }
            return new Image(NewImage);
        }
        #endregion
    }
}