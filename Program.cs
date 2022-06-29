//ClusterAlign main engine: Find and track fiducials in images of a tilt series
//Copyright (C) 2021-2022 by Shahar Seifer, Elabum lab, Weizmann Institute of Science
/*
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

using System;
using System.Text;
using System.IO;
using Emgu.CV;
using Emgu.Util;
using Emgu.CV.Util;
using Emgu.CV.Structure;
using Emgu.CV.CvEnum;
using System.Drawing;
using System.Threading;
using Emgu.CV.Features2D;
using System.Collections.Generic;

namespace ClusterAlign
{
    public class Program
    {
        static String Version_information = "ClusterAlign (ver 2022-June-29).";
        static bool xisRotation = ClusterAlign.Settings4ClusterAlign2.Default.xisRotation;
        static svector[] match_tolerance;
        static int cluster_size = ClusterAlign.Settings4ClusterAlign2.Default.cluster_size; //max radius of a single cluster
        static bool coswindow = ClusterAlign.Settings4ClusterAlign2.Default.coswindow; //coswindow: false- simple data, true- window already scalled by cos(tiltangle)= Hoppe sampling
        static double tolfactor_center = 0.01*ClusterAlign.Settings4ClusterAlign2.Default.TolFidCenter;// X fidsize= error in center location of fiducial
        static float fidsize_ratiovar = 0.01F*ClusterAlign.Settings4ClusterAlign2.Default.TolFidSize;  //ratio of variance to total size of fiducials
        static double[] z_values;
        static double[] tiltangles;
        static svector[,] svectors;
        static int[,] svectors_radius;
        static bool ignore_collisions = false;
        static bool cluster_visualize = false; //(Shows example of cluster visually at the end of the first iteration, depending on selected n1n1visualize)
        static int n1visualize = 336; ///336 -> check in advance a suitable numbers from fidn[ncenter,:] 
        static int[,,] visualize = new int[100, cluster_visualize? 50000:1, 2];//only for cluster visualization
        public static bool isMRCfile = true;

        public static void MyMain() 
        {
            Console.WriteLine("");
            Console.WriteLine(Version_information);
            Console.WriteLine("Written by Shahar Seifer, Elbaum lab, Weizmann Institute of Science.");
            Console.WriteLine("Note: GNU General Public License");
            string path = ClusterAlign.Settings4ClusterAlign2.Default.Path;
            string slash = "";
            if (path.Length>1)
            {
                if (path.Substring(path.Length-1)!="/" && path.Substring(path.Length - 1) != "\\")
                {
                    slash = (path.Contains("\\") ? "\\" : "/"); //depends on windows or linux
                }
            }
            string FileName = path + slash + ClusterAlign.Settings4ClusterAlign2.Default.DataFileName; //dataset, may be tif or mrc
            if (Path.GetExtension(FileName) == ".tif")
            {
                isMRCfile = false;
            }

            if (ClusterAlign.Settings4ClusterAlign2.Default.isArbAngle)
            {
                if (isMRCfile)
                {
                    string oldFileName = FileName;
                    FileName = path + slash + Path.GetFileNameWithoutExtension(FileName) + "_R.mrc";
                    Console.WriteLine("");
                    Console.WriteLine("Saving image with rotation axis parallel to Y (angle=0).");
                    double[] Dx_vect_zero = new double[360];
                    double[] Dy_vect_zero = new double[360];
                    double declared_phi = ClusterAlign.Settings4ClusterAlign2.Default.ArbAngle*Math.PI/180;
                    MrcStack myMrcStack = new MrcStack();
                    //will make the rotation reduced to zero, the new file will have rotation axis Y (phi=0), xisrotation is already false
                    myMrcStack.Fexport(oldFileName, FileName, ref Dx_vect_zero, ref Dy_vect_zero, ref declared_phi);

                }
                else
                {
                    Console.WriteLine("You request requires MRC file");
                    Thread.Sleep(10000);
                    System.Environment.Exit(0);
                }
            }
            string FidFileName = path + slash + Path.GetFileNameWithoutExtension(FileName)+ ".fid.txt";
            string NogapsFidFileName = path + slash + Path.GetFileNameWithoutExtension(FileName) + ".nogaps.fid.txt";
            string out_filename = path + slash + Path.GetFileNameWithoutExtension(FileName) + "_jali.mrc"; //our own ali files, just aligned mrc files, preferable if you use our reconstruction
            string normout_filename = path + slash + Path.GetFileNameWithoutExtension(FileName) + "_ali.mrc"; //both aligned and rotation axis in vertical direciton according to IMOD convention
            string showcluster_filename = path + slash + Path.GetFileNameWithoutExtension(FileName) + ".clusters.mrc";
            string report_filename = path + slash + Path.GetFileNameWithoutExtension(FileName) + ".output.txt";
            string basefilename = path + slash + Path.GetFileNameWithoutExtension(FileName);
            bool auto_fid_size = false;

            float fidsize = ClusterAlign.Settings4ClusterAlign2.Default.fidsize; //approx. fiducial size in pixels
            if (fidsize<=0)
            {
                auto_fid_size = true; //fidsize and masks will be updated before second iteration
                fidsize = 7; //tentative number for first iteration
            }
            int NfidMax = ClusterAlign.Settings4ClusterAlign2.Default.NfidMax; //maximum number of fiducials acquired in frame
            int N_minimum_tracked_fiducials_percent = ClusterAlign.Settings4ClusterAlign2.Default.N_minimum_tracked_fiducials; //percent of tracking useful for holding a fiducial
            double CorrectAngle = 0;
            double PreAlignmentTolx=ClusterAlign.Settings4ClusterAlign2.Default.PreAlignmentTol*(ClusterAlign.Settings4ClusterAlign2.Default.xisRotation? 1:1);
            double PreAlignmentToly = ClusterAlign.Settings4ClusterAlign2.Default.PreAlignmentTol * (ClusterAlign.Settings4ClusterAlign2.Default.xisRotation ? 1 : 1);
            int ncenter = ClusterAlign.Settings4ClusterAlign2.Default.ncenter;// best slice to see fiducials. if -1 then choose later:Nslices / 2
            float ThresholdG;//Threshold of "radial divergence" in gray map expected of fiducial
            int HKerSize = Convert.ToInt32(MathF.Ceiling(0.5F * fidsize * (1 + fidsize_ratiovar))); //Half kernel size for convolution in image processing
            int KerSize = 2 * HKerSize + 1;
            int HKerSizeDilate = (int)(HKerSize / 2);
            int KerSizeDilate = 2 * HKerSizeDilate + 1;
            int SptHKerSize = Convert.ToInt32(MathF.Ceiling(0.5F * fidsize )); //Half kernel size for convolution in image processing
            int Blursize = (SptHKerSize%2)==0? SptHKerSize+1: SptHKerSize;
            int SptKerSize = 2 * SptHKerSize + 1;
            int SptHKerSizeOut = (int)(SptHKerSize*MathF.Sqrt(2));
            int SptKerSizeOut = 2 * SptHKerSizeOut + 1;
            int HKerSize2nd = (int)(HKerSize*2); 
            int KerSize2nd = 2 * HKerSize2nd + 1;
            int HKerSize2ndOut = (int)(HKerSize * MathF.Sqrt(8)); 
            int KerSize2ndOut = 2 * HKerSize2ndOut + 1;
            double minVal = 0;
            double maxVal = 0;
            int exclude_radius = Math.Max(12, (int)(fidsize*1.2));
            //int level_expected = (int)(0.3*NfidMax); //number of cluster members expected, very tentative number
            System.Drawing.Point locationh = new System.Drawing.Point(0, 0);
            System.Drawing.Point locationl = new System.Drawing.Point(0, 0);
            bool fiducials_bright = ClusterAlign.Settings4ClusterAlign2.Default.fiducials_bright; //depends on the image B/D field. Fiducials are dark in bright mode so write false.
            const bool divergence_bright = true; //always true
            //int seedof3 =1 ; //if =1: seed of match is 3 vector match, otherwise 2 vector match.
            string[] tiltlines=null;
            bool isAttentionRequest = false;
            bool PseudoAttentionRequest = false;
            bool isfidRequest = false;
            int attention_NfidMax = Window1.AttentionNumMarkedfids;//Maximum number of fiducials selected by paintbrush in attention stack
            int attention_size = 0;
            NfidMax = NfidMax + attention_NfidMax;
            int NumofIterations = 1;
            int IterationNum;
            double[,] Bfinal;
            double[] Dx_vect;
            double[] Dy_vect;
            double report_phi = 0;
            double report_psi = 0;


            System.IO.FileStream reportoutfile = File.Create(report_filename);
            var reportfobj = new StreamWriter(reportoutfile, System.Text.Encoding.UTF8);
            reportfobj.WriteLine(Version_information);
            reportfobj.WriteLine("Tilt-series file: " + FileName);

            string rawtltFilename = @path + slash + ClusterAlign.Settings4ClusterAlign2.Default.TiltFileName;
            try
            {
               tiltlines = File.ReadAllLines(@path + slash + ClusterAlign.Settings4ClusterAlign2.Default.TiltFileName);  // datasetname.rawtlt
            }
            catch
            {
                Console.WriteLine("## The tilt file could not be found, check file name. ##");
                Thread.Sleep(10000);
                System.Environment.Exit(0);
            }
            tiltangles = new double[tiltlines.Length];
            for (int n = 0; n < tiltlines.Length; n++)
            {
                tiltangles[n] = (Convert.ToDouble(tiltlines[n])+CorrectAngle) * Math.PI / 180;
            }
            Mat[] slices = null;
            int Nslices = 0;
            int nslice;
            if (Path.GetExtension(FileName) == ".tif")
            {
                slices = CvInvoke.Imreadmulti(FileName, ImreadModes.Grayscale);
                Nslices = slices.Length;
                isMRCfile = false;
            }
            else if (Path.GetExtension(FileName) == ".mrc" || Path.GetExtension(FileName) == ".st" || Path.GetExtension(FileName) == ".preali")
            {
                MrcStack myMrcStack = new MrcStack();
                slices = myMrcStack.FOpen(FileName);
                Nslices = slices.Length;
                isMRCfile = true;
            }
            else
            {
                Console.WriteLine("File type name must be either mrc / tif / st / preali");
                Thread.Sleep(10000);
                System.Environment.Exit(0);
            }

            if (Nslices!= tiltlines.Length)
            {
                if (Nslices> tiltlines.Length)
                {
                    Nslices = tiltlines.Length;
                    Console.WriteLine("## Tilt-angle file not full, using part of the tilt-series. ##");
                }
                else
                {
                    Console.WriteLine("Too many lines in the tilt-angle file, using part of the list.");
                }
            }

            int[,] fidx = new int[Nslices, NfidMax]; 
            int[,] fidy = new int[Nslices, NfidMax]; 
            int[,] fidn = new int[Nslices, NfidMax]; 
            for (nslice=0; nslice<Nslices; nslice++ )
            {
                for (int n=0; n<NfidMax; n++)
                { fidn[nslice, n] = -1; }
            }
            int fid_count = 0;
            int[] match_score=new int[NfidMax];
            match_tolerance = new svector[Nslices];

            Mat[] Attentionslices=new Mat[Nslices];
            Mat[] PsAttentionslices = new Mat[Nslices];
            if (ClusterAlign.Window1.loadAttentionFile!="")
            { try
                {
                    Attentionslices = CvInvoke.Imreadmulti(ClusterAlign.Window1.loadAttentionFile, ImreadModes.Grayscale); //load tif stack similar in size to FileName
                    if (Attentionslices.Length!=Nslices)
                    {
                        Console.WriteLine("Attention file ignored due to incorrect format.");
                        ClusterAlign.Window1.loadAttentionFile = "";
                    }
                    else
                    {
                        isAttentionRequest = true;
                        reportfobj.WriteLine("Using attention file.");

                     }
                }
                catch
                {
                    Console.WriteLine("Attention file ignored due to incorrect format.");
                    ClusterAlign.Window1.loadAttentionFile = "";
                }
            }
            if (ClusterAlign.Window1.loadfidFile != "")
            {
                try
                {
                    readIMODfidModel(ref fidx, ref fidy, ref fidn, ref fid_count, ClusterAlign.Window1.loadfidFile);
                    isfidRequest = true;
                    attention_size = 10;
                    reportfobj.WriteLine("Imported fiducial file: "+ ClusterAlign.Window1.loadfidFile);
                }
                catch
                {
                    Console.WriteLine("fid file ignored due to incorrect format.");
                    ClusterAlign.Window1.loadfidFile = "";
                }
            }

            if (ncenter==-1)
            {
                ncenter = Nslices / 2;
                double mintilt = 1;
                for (int n = 0; n < tiltlines.Length; n++)
                {
                    if (Math.Abs(tiltangles[n])<mintilt)
                    {
                        mintilt = Math.Abs(tiltangles[n]);
                        ncenter = n;
                    }
                }

            }
            tp[,] locations = new tp[Nslices, NfidMax];
            int[] NFid = new int[Nslices];
            int Nrows = slices[0].Rows;
            int Ncols = slices[0].Cols;
            if (isMRCfile) //is going to be transposed
            {
                int Nstore= Nrows;
                Nrows = Ncols;
                Ncols = Nstore;
            }
            Mat slice_mat = new Mat(Nrows, Ncols, DepthType.Cv32F, 1);
            Mat Attention_slice_mat = new Mat(Nrows, Ncols, DepthType.Cv32F, 1);
            Mat matbuffer = new Mat(Nrows, Ncols, DepthType.Cv32F, 1);
            Mat matbuffer2 = new Mat(Nrows, Ncols, DepthType.Cv32F, 1);
            Mat matbuffer0 = new Mat(Nrows, Ncols, DepthType.Cv32F, 1);
            Mat zeros0 = new Mat(Nrows, Ncols, DepthType.Cv32F, 1);
            Mat ones0 = new Mat(Nrows, Ncols, DepthType.Cv32F, 1);
            ones0.SetTo(new MCvScalar(1.0));
            Mat localmaxima = new Mat(Nrows, Ncols, DepthType.Cv8U, 1);
            Mat matbuffer3 = new Mat(Nrows, Ncols, DepthType.Cv8U, 1);
            Mat matbuffer4 = new Mat(Nrows, Ncols, DepthType.Cv8U, 1);
            Mat matbuffer5 = new Mat(Nrows, Ncols, DepthType.Cv16U, 1);
            Mat matbuffer6 = new Mat(Nrows, Ncols, DepthType.Cv8U, 1);
            Mat matbuffer7 = new Mat(Nrows, Ncols, DepthType.Cv32F, 1);
            Mat matbuffer8 = new Mat(Nrows, Ncols, DepthType.Cv32F, 1);
            Array arr_matbuffer4 = new Byte[Nrows, Ncols];
            Array arr_matbuffer3 = new Byte[Nrows, Ncols];
            Mat gradx = new Mat(Nrows, Ncols, DepthType.Cv32F, 1);
            Mat grady = new Mat(Nrows, Ncols, DepthType.Cv32F, 1);
            int kmask_area = 0;
            float kradius;
            Mat submatbuffer; //= new Mat(9, 9, DepthType.Cv32F, 1);
            Mat submatbuffer_match = new Mat(Nrows-8, Ncols-8, DepthType.Cv32F, 1);
            Mat submatbuffer_match_attention = new Mat(Nrows - 8, Ncols - 8, DepthType.Cv32F, 1);
            Mat grand_submatbuffer= new Mat(exclude_radius * 2 + 1, exclude_radius * 2 + 1, DepthType.Cv32F, 1); ;
            grand_submatbuffer.SetTo(new MCvScalar(0));

            IntPtr kerx_ptr = new IntPtr();
            IntPtr kery_ptr = new IntPtr();
            IntPtr kerasymx_ptr = new IntPtr();
            IntPtr kerasymy_ptr = new IntPtr();
            IntPtr kerdiag1_ptr = new IntPtr();
            IntPtr kerdiag2_ptr = new IntPtr();
            IntPtr keradial_ptr = new IntPtr();
            kerx_ptr = CvInvoke.cvCreateMat(KerSize, KerSize, DepthType.Cv32F);
            kery_ptr = CvInvoke.cvCreateMat(KerSize, KerSize, DepthType.Cv32F);
            kerasymx_ptr = CvInvoke.cvCreateMat(KerSize, KerSize, DepthType.Cv32F);
            kerasymy_ptr = CvInvoke.cvCreateMat(KerSize, KerSize, DepthType.Cv32F);
            kerdiag1_ptr = CvInvoke.cvCreateMat(KerSize, KerSize, DepthType.Cv32F);
            kerdiag2_ptr = CvInvoke.cvCreateMat(KerSize, KerSize, DepthType.Cv32F);
            keradial_ptr = CvInvoke.cvCreateMat(KerSize, KerSize, DepthType.Cv32F);
            for (int xind = -HKerSize; xind <= HKerSize; xind++)
            {
                for (int yind = -HKerSize; yind <= HKerSize; yind++)
                {
                    kradius = MathF.Sqrt(xind * xind + yind * yind) + 0.1F;
                    if (kradius >= 0.5F * fidsize * (1 - fidsize_ratiovar) && kradius <= 0.5F * fidsize * (1 + fidsize_ratiovar))
                    {
                        CvInvoke.cvSetReal2D(kerx_ptr, yind + HKerSize, xind + HKerSize, (float)xind / kradius);
                        CvInvoke.cvSetReal2D(kery_ptr, yind + HKerSize, xind + HKerSize, (float)yind / kradius);
                    }
                    else
                    {
                        CvInvoke.cvSetReal2D(kerx_ptr, yind + HKerSize, xind + HKerSize, 0);
                        CvInvoke.cvSetReal2D(kery_ptr, yind + HKerSize, xind + HKerSize, 0);
                    }
                    CvInvoke.cvSetReal2D(kerasymx_ptr, yind + HKerSize, xind + HKerSize, (float)MathF.Sign(xind));
                    CvInvoke.cvSetReal2D(kerasymy_ptr, yind + HKerSize, xind + HKerSize, (float)MathF.Sign(yind));
                    CvInvoke.cvSetReal2D(kerdiag1_ptr, yind + HKerSize, xind + HKerSize, (float)MathF.Sign(0.5f * xind + 0.5f * yind));
                    CvInvoke.cvSetReal2D(kerdiag2_ptr, yind + HKerSize, xind + HKerSize, (float)MathF.Sign(0.5f * xind - 0.5f * yind));
                    CvInvoke.cvSetReal2D(keradial_ptr, yind + HKerSize, xind + HKerSize, MathF.Sqrt((float)(MathF.Pow(xind,2)+MathF.Pow(yind,2))));
                }
            }
            Mat kerx = CvInvoke.CvArrToMat(kerx_ptr);
            Mat kery = CvInvoke.CvArrToMat(kery_ptr);
            Mat kerasymx = CvInvoke.CvArrToMat(kerasymx_ptr);
            Mat kerasymy = CvInvoke.CvArrToMat(kerasymy_ptr);
            Mat kerdiag1 = CvInvoke.CvArrToMat(kerdiag1_ptr);
            Mat kerdiag2 = CvInvoke.CvArrToMat(kerdiag2_ptr);
            Mat keradial = CvInvoke.CvArrToMat(keradial_ptr);
            Mat ex_kmask = CvInvoke.GetStructuringElement(ElementShape.Ellipse, new Size(KerSize, KerSize), new Point(HKerSize, HKerSize));
            Mat ex_kmaskDilate = CvInvoke.GetStructuringElement(ElementShape.Ellipse, new Size(KerSizeDilate, KerSizeDilate), new Point(HKerSizeDilate, HKerSizeDilate));
            Mat ex_kmask2nd = CvInvoke.GetStructuringElement(ElementShape.Ellipse, new Size(KerSize2nd, KerSize2nd), new Point(HKerSize2nd, HKerSize2nd));
            Mat ex_kmask2ndOut = CvInvoke.GetStructuringElement(ElementShape.Ellipse, new Size(KerSize2ndOut, KerSize2ndOut), new Point(HKerSize2ndOut, HKerSize2ndOut));
            Mat ex_Sptkmask = CvInvoke.GetStructuringElement(ElementShape.Ellipse, new Size(SptKerSize, SptKerSize), new Point(SptHKerSize, SptHKerSize));
            Mat ex_SptkmaskOut = CvInvoke.GetStructuringElement(ElementShape.Ellipse, new Size(SptKerSizeOut, SptKerSizeOut), new Point(SptHKerSizeOut, SptHKerSizeOut));
            Mat ex_symask = CvInvoke.GetStructuringElement(ElementShape.Ellipse, new Size(KerSize, KerSize), new Point(HKerSize, HKerSize));
            Mat ex_3by3 = CvInvoke.GetStructuringElement(ElementShape.Ellipse, new Size(3, 3), new Point(2, 2));
            Mat ex_kerx = new Mat(KerSize, KerSize, DepthType.Cv32F, 1); 
            Mat ex_kery = new Mat(KerSize, KerSize, DepthType.Cv32F, 1);
            Size sayzero = new Size(0, 0);
            kmask_area =(int)(3.14F* HKerSize* HKerSize);
            int ex_kmask_area;
            int high_pass_size;
            string win1;
            svectors = new svector[Nslices, NfidMax * NfidMax];
            svectors_radius = new int[Nslices, NfidMax * NfidMax];
            float tolerance_radius = 0;
            int[] Nsvector = new int[Nslices];
            float col1, col2, row1, row2;
            float dx_norm, dy_norm;
            svector svecnull;
            svecnull = new svector(0,0,0,0,0,0);
            int NumberofLevels = ClusterAlign.Settings4ClusterAlign2.Default.Ncluster - 1; //number of neccessary matching vectors between cluster instances (Ncluster number of necessary fiducials in cluster)
            int[] cluster_candidate_list = new int[NfidMax];
            int[,] sfid_cluster_list = new int[Nslices, NfidMax];  //list of elements numbers in each slice participating in the current gathered cluster
            int[,] smatch = new int[Nslices, NfidMax]; // correspondance of fiducial numbers in each slice to the fiducial in center slice
            int[,] smatch_strength1 = new int[Nslices, NfidMax]; // strength of matching
            int[,] smatch_strength2 = new int[Nslices, NfidMax]; // strength of matching
            int[] nextsv = new int[Nslices];
            int[] match_counter = new int[NfidMax];
            int minslice=0;
            int maxslice=0;
            int count;
            bool missing_location;
            int lastpx = 0;
            int lastpy = 0;
            Dx_vect = new double[Nslices];
            Dy_vect = new double[Nslices];
            Array.Clear(Dx_vect, 0, Nslices);
            Array.Clear(Dy_vect, 0, Nslices);
            double xc_ncenter=Ncols / 2;
            double yc_ncenter=Nrows / 2;
            double xc= Ncols / 2;
            double yc= Nrows / 2;
            int xc_ns, yc_ns;
            int marginx = 0;
            int marginy = 0;
            int grand_acc_count = 0;
            int blursizex, blursizey;
            double noise_std = 32000;


            if (ClusterAlign.Settings4ClusterAlign2.Default.ForceFill)
            {
                NumofIterations = 2;
            }
            bool optical_test = true;
  

            double maxz = 0.05 * Nrows; //maximum z differences we might handle with this software
            for (nslice = 0; nslice < Nslices; nslice++)
            {
                match_tolerance[nslice] = new svector(0, 0, (int)Math.Floor(tolfactor_center * fidsize), (int)Math.Floor(tolfactor_center * fidsize), (int)(Math.Max(PreAlignmentTolx,150) * 2 + maxz * Math.Sin(Math.Abs(tiltangles[nslice]))), (int)(Math.Max(PreAlignmentToly,150) * 2 + maxz * Math.Sin(Math.Abs(tiltangles[nslice]))));
            }
            
            for (IterationNum = 0; IterationNum < NumofIterations; IterationNum++)
            {

                if (isfidRequest)
                {
                    if (IterationNum == 0)
                        optical_test = false;//skip optical recognition, just to use loaded fid points
                    else
                    {
                        optical_test = true;
                        isfidRequest = false;
                    }
                }

                tolerance_radius = (float)(tolfactor_center * fidsize);
                win1 = "Fiducials Marked";
                CvInvoke.NamedWindow(win1, WindowFlags.KeepRatio);
                for (nslice = 0; nslice < Nslices; nslice++)
                {

                    slices[nslice].ConvertTo(slice_mat, DepthType.Cv32F, 1, 0);
                    if (isMRCfile)
                    {
                        CvInvoke.Transpose(slice_mat, slice_mat); //switch x and y axis so afterwards x(col) is along 3dmod x and y[row] is along 3dmod y. The origin is upperleft corner here and is lowerleft corner in 3dmod but it is the same relative to the image
                    }
                    else
                    {
                        CvInvoke.Flip(slice_mat, slice_mat, FlipType.Vertical);
                    }

                    MCvScalar find_mean = new MCvScalar(0);
                    MCvScalar find_std = new MCvScalar(0);

                    // Now x/col axis is to the right and y/row axis is downward, starting at uperleftcorner in this program. Up and down are reversed compared to 3dmod but the y values are consistent with the image
                    int lp_size = (int)Math.Round(fidsize / 3.0);
                    if (lp_size % 2 == 0) lp_size--;
                    if (lp_size < 1) lp_size = 1;
                    //Low-pass filter
                    CvInvoke.GaussianBlur(slice_mat, slice_mat, new System.Drawing.Size(lp_size, lp_size), 0, 0);

                    if (optical_test)
                    {

                        //The size of spherical fiducials will expand in one direction in the consine sampling aspect ratio
                        //Note that the coordinates to Mat array go like y,x but in some commands the order is x,y 
                        double say1 = 1d /(Math.Cos(tiltangles[nslice]));
                        double say1dx = (!coswindow) ? 1d : (xisRotation? 1d: say1);
                        double say1dy = (!coswindow) ? 1d : (xisRotation ? say1 : 1d); 
                        CvInvoke.Resize(kerx, ex_kerx, sayzero, say1dx, say1dy, Inter.Cubic); //resize the mat file ex_kerx to say1dy,say1dx 
                        CvInvoke.Resize(kery, ex_kery, sayzero, say1dx, say1dy, Inter.Cubic);
                        ex_kmask = CvInvoke.GetStructuringElement(ElementShape.Ellipse, new Size((int)(KerSize * say1dx) + (int)(KerSize * say1dx + 1) % 2, (int)(KerSize * say1dy) + (int)(KerSize * say1dy + 1) % 2), new Point(-1, -1));
                        ex_kmask_area = CvInvoke.CountNonZero(ex_kmask);
                        ex_kmaskDilate = CvInvoke.GetStructuringElement(ElementShape.Ellipse, new Size((int)(KerSizeDilate * say1dx) + (int)(KerSizeDilate * say1dx + 1) % 2, (int)(KerSizeDilate * say1dy) + (int)(KerSizeDilate * say1dy + 1) % 2), new Point(-1, -1));

                        ex_kmask2nd = CvInvoke.GetStructuringElement(ElementShape.Ellipse, new Size((int)(KerSize2nd * say1dx) + (int)(KerSize2nd * say1dx + 1) % 2, (int)(KerSize2nd * say1dy) + (int)(KerSize2nd * say1dy + 1) % 2), new Point(-1, -1));
                        ex_kmask2ndOut = CvInvoke.GetStructuringElement(ElementShape.Ellipse, new Size((int)(KerSize2ndOut * say1dx) + (int)(KerSize2ndOut * say1dx + 1) % 2, (int)(KerSize2ndOut * say1dy) + (int)(KerSize2ndOut * say1dy + 1) % 2), new Point(-1, -1));
                        ex_Sptkmask = CvInvoke.GetStructuringElement(ElementShape.Ellipse, new Size((int)(SptKerSize * say1dx) + (int)(SptKerSize * say1dx + 1) % 2, (int)(SptKerSize * say1dy) + (int)(SptKerSize * say1dy + 1) % 2), new Point(-1, -1));
                        ex_SptkmaskOut = CvInvoke.GetStructuringElement(ElementShape.Ellipse, new Size((int)(SptKerSizeOut * say1dx) + (int)(SptKerSizeOut * say1dx + 1) % 2, (int)(SptKerSizeOut * say1dy) + (int)(SptKerSizeOut * say1dy + 1) % 2), new Point(-1, -1));


                        if (nslice == ncenter)
                        {
                            CvInvoke.MeanStdDev(slice_mat, ref find_mean, ref find_std);
                            noise_std = find_std.V0;
                        }

                        //High pass filter to matbuffer
                        high_pass_size = (int)MathF.Floor(fidsize*20 );
                        high_pass_size = (high_pass_size % 2 == 1) ? high_pass_size : high_pass_size + 1;
                        CvInvoke.GaussianBlur(slice_mat, matbuffer2, new System.Drawing.Size(high_pass_size, high_pass_size), 0, 0);
                        CvInvoke.Subtract(slice_mat, matbuffer2, matbuffer); //remove large features (like shadows), store in matbuffer
                                                                             
                        blursizex = (int)MathF.Floor((float)say1dx * fidsize/2 );
                        blursizey = (int)MathF.Floor((float)say1dy * fidsize/2 );
                        blursizex = blursizex % 2 == 1 ? blursizex : blursizex + 1; //must be odd numbers
                        blursizey = blursizey % 2 == 1 ? blursizey : blursizey + 1;

                        //win1 = "normal";
                        //Program.show_grayimage(slice_mat, win1, Nrows, Ncols);

                        //**** acquire gradients in the x and y directions ****
                        Program.SpatialGradient(slice_mat, ref gradx, ref grady); //must use non high pass filtered so the margins will make 0 divergence
                        // **** Calculate Sqrt((Gradx(*)Kerx)^2+(Grady(*)Kery)^2).  (*) denotes convolution
                        CvInvoke.Filter2D(gradx, matbuffer0, ex_kerx, anchor: new System.Drawing.Point(-1, -1));// Convolute gradx with kerx,  anchor point at the center , save in matbuffer0
                        CvInvoke.Multiply(matbuffer0, matbuffer0, matbuffer2);//matbuffer2=matbuffer0.^2
                        CvInvoke.Filter2D(grady, matbuffer0, ex_kery, anchor: new System.Drawing.Point(-1, -1));// Convolute grady with kery,  anchor point at the center , save in matbuffer0
                        CvInvoke.Multiply(matbuffer0, matbuffer0, matbuffer0);//matbuffer0=matbuffer0.^2
                        CvInvoke.Add(matbuffer2, matbuffer0, matbuffer2);//matbuffer2=matbuffer2+matbuffer0
                        CvInvoke.Sqrt(matbuffer2, matbuffer7);//matbuffer7=sqrt(matbuffer2) = divergence
                        //remove high divergence lines due to borders
                        Image<Gray, float> img7 = matbuffer7.ToImage<Gray, float>();
                        for (int linn = 2; linn < Ncols / 3; linn++)
                        {
                            img7.ROI = new Rectangle(1, 1, linn - 1, Nrows - 1);
                            if (img7.CountNonzero()[0] > 0)
                            {
                                img7.ROI = new Rectangle(0, 0, linn + 12, Nrows);
                                img7.SetValue(new MCvScalar(0));
                                break;
                            }
                        }
                        for (int linn = 2; linn < Nrows / 3; linn++)
                        {
                            img7.ROI = new Rectangle(1, 1, Ncols - 2, linn - 1);
                            int cntr = img7.CountNonzero()[0];
                            if (img7.CountNonzero()[0] > 0)
                            {
                                img7.ROI = new Rectangle(0, 0, Ncols, linn + 12);
                                img7.SetValue(new MCvScalar(0));
                                break;
                            }
                        }
                        for (int linn = 2; linn < Ncols / 3; linn++)
                        {
                            img7.ROI = new Rectangle(Ncols - linn, 1, linn - 1, Nrows - 2);
                            if (img7.CountNonzero()[0] > 0)
                            {
                                img7.ROI = new Rectangle(Ncols - linn - 12, 0, linn + 12, Nrows);
                                img7.SetValue(new MCvScalar(0));
                               break;
                            }
                        }
                        for (int linn = 2; linn < Nrows / 3; linn++)
                        {
                            img7.ROI = new Rectangle(1, Nrows - linn, Ncols - 2, linn - 1);
                            if (img7.CountNonzero()[0] > 0)
                            {
                                img7.ROI = new Rectangle(0, Nrows - linn - 12, Ncols, linn + 12);
                                img7.SetValue(new MCvScalar(0));
                                break;
                            }
                        }
                        img7.ROI = Rectangle.Empty;
                        matbuffer7 = img7.Mat;   //This is image sensitive to gradients around the markers that fit inside the ring of expected radii
                        CvInvoke.GaussianBlur(matbuffer7, matbuffer7, new System.Drawing.Size(3, 3), 0, 0);
                        //win1 = "mat7=div";
                        //Program.show_grayimage(matbuffer7, win1, Nrows, Ncols);
                        //CvInvoke.WaitKey(2000);

                        if (isAttentionRequest)
                        {
                            //painted with color 0 over important fiducials - > will generate mask in matbuffer6
                            Attentionslices[nslice].ConvertTo(Attention_slice_mat, DepthType.Cv8U, 1, 0);
                            CvInvoke.Flip(Attention_slice_mat, Attention_slice_mat, FlipType.Vertical);
                            CvInvoke.Threshold(Attention_slice_mat, matbuffer6, 1, 1, ThresholdType.BinaryInv); //matbuffer6= Attention_slice_mat<1? 1:0 (used as mask to pick color 0)
                            CvInvoke.Multiply(matbuffer7, matbuffer6, matbuffer7, 1.0, DepthType.Cv32F);
                            CvInvoke.Multiply(matbuffer, matbuffer6, matbuffer, 1.0, DepthType.Cv32F);
                        }
                        //keep matbuffer6 that holds attention mask, and update matbuffer7 to include features only at the attention spots
                        // win1 = "mat7=attention";
                        // Program.show_grayimage(matbuffer7, win1, Nrows, Ncols);
                        //win1 = "mat=attention";
                        //Program.show_grayimage(matbuffer, win1, Nrows, Ncols);




                        CvInvoke.MinMaxLoc(matbuffer7, ref minVal, ref maxVal, ref locationl, ref locationh);
                        matbuffer7.ConvertTo(matbuffer7, DepthType.Cv32F, 1d, -minVal);
                        
                        if (!fiducials_bright)
                        { CvInvoke.Subtract(zeros0, matbuffer, matbuffer); }
                        CvInvoke.MinMaxLoc(matbuffer, ref minVal, ref maxVal, ref locationl, ref locationh);
                        matbuffer.ConvertTo(matbuffer, DepthType.Cv32F, 1d, -minVal);
                        CvInvoke.GaussianBlur(matbuffer, matbuffer, new System.Drawing.Size(3, 3), 0, 0);
                        CvInvoke.Multiply(matbuffer, matbuffer, matbuffer0);//matbuffer0=matbuffer*matbuffer

                        MCvScalar avgmat = CvInvoke.Mean(matbuffer0, null);
                        if (!ClusterAlign.Settings4ClusterAlign2.Default.coswindow)
                        {
                            if (ClusterAlign.Settings4ClusterAlign2.Default.xisRotation)
                            {
                                int ysize = (int)((1.0 - MathF.Cos((float)tiltangles[nslice])) * Nrows / 2.0);
                                CvInvoke.Rectangle(matbuffer0, new Rectangle(0, 0, Ncols, ysize), avgmat, -1); //(x,y,width,height), upperleft corner, color=0, linewidth -1: filled rectangle
                                CvInvoke.Rectangle(matbuffer0, new Rectangle(0, Nrows - ysize, Ncols, ysize), avgmat, -1); //(x,y,width,height), upperleft corner, color=0, linewidth -1: filled rectangle
                                marginy = ysize;
                                marginx = 0;
                            }
                            else
                            {
                                int xsize = (int)((1.0 - MathF.Cos((float)tiltangles[nslice])) * Ncols / 2.0);
                                CvInvoke.Rectangle(matbuffer0, new Rectangle(0, 0, xsize, Nrows), avgmat, -1); //(x,y,width,height), upperleft corner, color=0, linewidth -1: filled rectangle
                                CvInvoke.Rectangle(matbuffer0, new Rectangle(Ncols - xsize, 0, xsize, Nrows), avgmat, -1); //(x,y,width,height), upperleft corner, color=0, linewidth -1: filled rectangle
                                marginy = 0;
                                marginx = xsize;
                            }
                        }
                        //show win1 = "mat*mat7";
                        //show Program.show_grayimage(matbuffer0, win1, Nrows, Ncols);

                        CvInvoke.Dilate(matbuffer7, matbuffer2, ex_kmaskDilate, anchor: new System.Drawing.Point(-1, -1), 1, BorderType.Default, new MCvScalar(0));//Dilate one iteration:  replace neighborhoods with maxima and save in matbuffer2
                        CvInvoke.Compare(matbuffer2, matbuffer7, localmaxima, CmpType.Equal); //localmaxima= xFF at locations of local peaks of divergenece (based on comparison of dilation with original image)
                        //win1 = "mat2";
                        //Program.show_grayimage(matbuffer2, win1, Nrows, Ncols);
                        //remove localmaxima that are actually zeros
                        CvInvoke.Compare(matbuffer2, zeros0, matbuffer3, CmpType.Equal);
                        CvInvoke.BitwiseXor(localmaxima, matbuffer3, localmaxima);

                        localmaxima.ConvertTo(matbuffer7, DepthType.Cv32F,1d/255);
                        CvInvoke.Dilate(matbuffer7, matbuffer7, ex_kmask, anchor: new System.Drawing.Point(-1, -1), 1, BorderType.Default, new MCvScalar(0));

                        if (isAttentionRequest)
                        {
                            //Use matbuffer6
                        }
                        else
                        {
                            matbuffer7.ConvertTo(matbuffer6, DepthType.Cv8U, 1, 0);
                        }
                        ThresholdG = Program.FindThreshold(matbuffer0, matbuffer6, NfidMax, ex_kmask_area, divergence_bright); //Find ThresoldG according to histogram of dilated divergence map
                        CvInvoke.Threshold(matbuffer0, matbuffer3, ThresholdG, 0xFF, ThresholdType.Binary); //matbuffer3= matbuffer4>ThresholdG? 0xFF:0
                        matbuffer3.ConvertTo(matbuffer3, DepthType.Cv8U);
                        if (fidsize>=15)    CvInvoke.Erode(matbuffer3, matbuffer3, ex_3by3, anchor: new System.Drawing.Point(-1, -1),(int)Math.Round(fidsize/15), BorderType.Default, new MCvScalar(0));
                        CvInvoke.BitwiseAnd(localmaxima, matbuffer3, matbuffer4); //matbuffer4= localmaxima & matbuffer3

                         bool flag_record;
                        {
                            Mat ROItemp;
                            Mat t_ROItemp = new Mat(exclude_radius * 2 + 1, exclude_radius * 2 + 1, DepthType.Cv32F, 1);
                            Mat scoreImg = new Mat(1, 1, DepthType.Cv32F, 1);
                            flag_record = false;
                            submatbuffer = new Mat(exclude_radius * 2 + 1, exclude_radius * 2 + 1, DepthType.Cv32F, 1);
                            submatbuffer.SetTo(new MCvScalar(0));
                            NFid[nslice] = 0;
                            arr_matbuffer4 = matbuffer4.GetData();  //from mat to array, stores the fiducials according to general procedure
                            for (int row = marginy + exclude_radius + 1; row < Nrows - marginy - exclude_radius; row++)
                            {
                                for (int col = marginx + exclude_radius + 1; col < Ncols - marginx - exclude_radius; col++)
                                {
                                    if (Convert.ToByte(arr_matbuffer4.GetValue(row, col)) == 0xFF && NFid[nslice] < NfidMax)
                                    {
                                        flag_record = true;
                                        for (int ind_exc = NFid[nslice] - 1; ind_exc >= 0; ind_exc--)
                                        {
                                            if (Math.Abs(row - locations[nslice, ind_exc].row) <= exclude_radius && Math.Abs(col - locations[nslice, ind_exc].col) <= exclude_radius)
                                            {
                                                if (Math.Sqrt(Math.Pow((row - locations[nslice, ind_exc].row), 2) + Math.Pow((col - locations[nslice, ind_exc].col), 2)) <= exclude_radius)
                                                {
                                                    flag_record = false;
                                                    break;
                                                }
                                            }
                                        }
                                        if (flag_record)
                                        {
                                            ROItemp = new Mat(matbuffer0, new Rectangle(col - exclude_radius, row - exclude_radius, exclude_radius * 2 + 1, exclude_radius * 2 + 1));
                                            //do not accumulate very anostropic features
                                            CvInvoke.Transpose(ROItemp, t_ROItemp);
                                            CvInvoke.MatchTemplate(ROItemp, t_ROItemp, scoreImg, TemplateMatchingType.CcorrNormed);
                                            maxVal = (double)scoreImg.GetValue(0, 0);
                                            if (maxVal > 0.4)
                                            {
                                                CvInvoke.Accumulate(ROItemp, submatbuffer);//find average subimage of the fiducials
                                                locations[nslice, NFid[nslice]] = new tp(row, col);// write in table the location of fiducial found (y,x)
                                                NFid[nslice]++;
                                            }
                                            else
                                            {
                                                //black out anisotropic locations
                                                CvInvoke.Circle(matbuffer0, new Point(col, row), exclude_radius, new MCvScalar(0), -1);
                                            }
                                        }
                                    }
                                }
                            }
                            if (NFid[nslice] > 20) 
                            { submatbuffer.ConvertTo(submatbuffer, DepthType.Cv32F, 1.0d / NFid[nslice]); } //normalize accumulated images to average 8U/32F
                            else
                            {
                                //intial assumption about shape of fiducial (artifical template):
                                submatbuffer.SetTo(new MCvScalar(0));
                                CvInvoke.Circle(submatbuffer, new Point(exclude_radius, exclude_radius), (int)(0.75*HKerSize), avgmat, -1); //avgmat instead of new MCvScalar();
                                CvInvoke.GaussianBlur(submatbuffer, submatbuffer, new System.Drawing.Size(blursizex, blursizey), 0, 0);
                            }
                            win1 = "generated template";
                            Program.show_grayimage(submatbuffer, win1, exclude_radius * 2 + 2, exclude_radius * 2 + 2);
                            if (IterationNum == 0 && Math.Abs(tiltangles[nslice]) < 30 * Math.PI / 180)
                            {
                                CvInvoke.Accumulate(submatbuffer, grand_submatbuffer);
                                grand_acc_count++;
                            }
                            if (IterationNum == 0 && nslice == Nslices - 1 && grand_acc_count>0) //operations in the last slice
                            {
                                grand_submatbuffer.ConvertTo(grand_submatbuffer, DepthType.Cv32F, 1.0d / grand_acc_count);//normalize
                                //CvInvoke.MeanStdDev(grand_submatbuffer,ref find_mean, ref find_std);
                                CvInvoke.MinMaxLoc(grand_submatbuffer, ref minVal, ref maxVal, ref locationl, ref locationh);
                                double signal_var = Math.Sqrt(maxVal)- Math.Sqrt(minVal); //this is the isolated signal (fiducial) with respect to the non-squared pixel levels
                                Console.WriteLine("CNR= {0:0.000}", signal_var/noise_std); //Contrast to noise ratio of the fiducials
                                reportfobj.WriteLine("CNR= {0:0.000}", signal_var / noise_std);
                                if (auto_fid_size && maxVal - minVal > 0)
                                {
                                    //this block adjusts fidsize and its related mask 
                                    double threshold_t = 0;
                                    Mat thresh = new Mat(grand_submatbuffer.Rows, grand_submatbuffer.Cols, DepthType.Cv8U, 1);
                                    grand_submatbuffer.ConvertTo(grand_submatbuffer, DepthType.Cv32F, 255d / (maxVal - minVal), -minVal * 255d / (maxVal - minVal));
                                    grand_submatbuffer.ConvertTo(matbuffer3, DepthType.Cv8U);
                                    threshold_t = CvInvoke.Threshold(matbuffer3, matbuffer3, 0, 255, ThresholdType.Otsu); //Finds sharp edge representation of the average fiducial by Otsu algorithm based on graidents
                                    MCvScalar sumim = CvInvoke.Sum(thresh);
                                    //Update avergae fiducial radius for next iteration
                                    fidsize = (float)(2*Math.Sqrt((sumim.V0 / 255.0) / Math.PI)); //based on area of a circle
                                    Console.WriteLine("Updated avg. fiducial size= {0:0.0}",fidsize);
                                    HKerSize = Convert.ToInt32(MathF.Ceiling(0.5F * fidsize * (1 + fidsize_ratiovar))); //Half kernel size for convolution in image processing
                                    KerSize = 2 * HKerSize + 1;
                                    HKerSizeDilate = (int)(HKerSize / 2);
                                    KerSizeDilate = 2 * HKerSizeDilate + 1;
                                    kerx_ptr = CvInvoke.cvCreateMat(KerSize, KerSize, DepthType.Cv32F);
                                    kery_ptr = CvInvoke.cvCreateMat(KerSize, KerSize, DepthType.Cv32F);
                                    for (int xind = -HKerSize; xind <= HKerSize; xind++)
                                    {
                                        for (int yind = -HKerSize; yind <= HKerSize; yind++)
                                        {
                                            kradius = MathF.Sqrt(xind * xind + yind * yind) + 0.1F;
                                            if (kradius >= 0.5F * fidsize * (1 - fidsize_ratiovar) && kradius <= 0.5F * fidsize * (1 + fidsize_ratiovar))
                                            {
                                                CvInvoke.cvSetReal2D(kerx_ptr, yind + HKerSize, xind + HKerSize, (float)xind / kradius);
                                                CvInvoke.cvSetReal2D(kery_ptr, yind + HKerSize, xind + HKerSize, (float)yind / kradius);
                                            }
                                            else
                                            {
                                                CvInvoke.cvSetReal2D(kerx_ptr, yind + HKerSize, xind + HKerSize, 0);
                                                CvInvoke.cvSetReal2D(kery_ptr, yind + HKerSize, xind + HKerSize, 0);
                                            }
                                        }
                                    }
                                    //Update to have better screening for candidates of template averge based on the gradient method used to locate fiducial positions
                                    kerx = CvInvoke.CvArrToMat(kerx_ptr);
                                    kery = CvInvoke.CvArrToMat(kery_ptr);
                                }
                            }

                        } 
 
                        CvInvoke.MatchTemplate(matbuffer0, submatbuffer, submatbuffer_match, TemplateMatchingType.CcorrNormed);//CcoeffNormed correlation of signal above average
                        if (PseudoAttentionRequest)
                        {
                            //painted with color 0 over circles of expected fiducial locations according to fit
                            PsAttentionslices[nslice].ConvertTo(Attention_slice_mat, DepthType.Cv8U, 1, 0);
                            CvInvoke.Threshold(Attention_slice_mat, matbuffer4, 1, 0xFF, ThresholdType.BinaryInv); //matbuffer6= Attention_slice_mat<1? 0xFF:0 (used as mask to pick color 0)
                            CvInvoke.Multiply(matbuffer4, matbuffer0, matbuffer2, 1.0, DepthType.Cv32F);
                            CvInvoke.MatchTemplate(matbuffer2, submatbuffer, submatbuffer_match_attention, TemplateMatchingType.CcorrNormed);//CcoeffNormed correlation of signal above average, with matbuffer6 as mask
                            CvInvoke.Resize(matbuffer4, matbuffer4, new Size(submatbuffer_match_attention.Width, submatbuffer_match_attention.Height));
                            double background_mask = (double)CvInvoke.Mean(submatbuffer_match_attention, matbuffer4).V0;
                            submatbuffer_match_attention.ConvertTo(submatbuffer_match_attention, DepthType.Cv32F, 1.0d, -background_mask); //subtract background level
                            Mat zerosa = new Mat(submatbuffer_match_attention.Height, submatbuffer_match_attention.Width, DepthType.Cv32F, 1);
                            CvInvoke.Max(zerosa, submatbuffer_match_attention, submatbuffer_match_attention);
                            CvInvoke.Add(submatbuffer_match, submatbuffer_match_attention, submatbuffer_match);
                        }

                        CvInvoke.GaussianBlur(submatbuffer_match, submatbuffer_match, new System.Drawing.Size(3, 3), 0, 0);
                        ///Console.WriteLine("Slice#="+nslice.ToString()+"  Number of fiducials by threshold="+ NFid[nslice].ToString());
                        ///win1 = "matching template";
                        ///Program.show_grayimage(submatbuffer_match, win1, Nrows- exclude_radius*2, Ncols- exclude_radius*2);
                        matbuffer2.SetTo(new MCvScalar(0));
                        submatbuffer_match.CopyTo(new Mat(matbuffer2, new Rectangle(exclude_radius, exclude_radius, submatbuffer_match.Cols, submatbuffer_match.Rows))); //matbuffer0 is now the image of correlation with submat 

                        CvInvoke.Dilate(matbuffer2, matbuffer7, ex_kmaskDilate, anchor: new System.Drawing.Point(-1, -1), 1, BorderType.Default, new MCvScalar(0));//Dilate one iteration:  replace neighborhoods with maxima and save in matbuffer2
                        CvInvoke.Compare(matbuffer2, matbuffer7, localmaxima, CmpType.Equal); //localmaxima= xFF at locations of local peaks of divergenece (based on comparison of dilation with original image)
                        localmaxima.ConvertTo(matbuffer7, DepthType.Cv32F,1d/255);
                        CvInvoke.Dilate(matbuffer7, matbuffer7, ex_kmask, anchor: new System.Drawing.Point(-1, -1), 1, BorderType.Default, new MCvScalar(0));
                        //
                        CvInvoke.Multiply(matbuffer0, matbuffer0, matbuffer0);//matbuffer0=matbuffer0*matbuffer0   //NEEDED TO EMPHASIS MORE THE ACTUAL IMAGE (maybe use filter on match image)
                        CvInvoke.Multiply(matbuffer2, matbuffer0, matbuffer0);//matbuffer0=matbuffer2*matbuffer0
                        CvInvoke.GaussianBlur(matbuffer0, matbuffer0, new System.Drawing.Size(lp_size, lp_size), 0, 0);
                        //string win3 = "for final threshold";
                        //CvInvoke.NamedWindow(win3, WindowFlags.KeepRatio);
                        //Program.show_grayimage(matbuffer0, win3, Nrows, Ncols);

                        //Calculate improved matbuffer3 and use it to caluclate matbuffer4 based on previous found localmax
                        if (isAttentionRequest)
                        { 
                            //Use matbuffer6
                        } 
                        else
                        {
                            matbuffer7.ConvertTo(matbuffer6, DepthType.Cv8U, 1, 0);
                        }
                        ThresholdG = Program.FindThreshold(matbuffer0, matbuffer6, NfidMax, ex_kmask_area, divergence_bright);
                        CvInvoke.Threshold(matbuffer0, matbuffer2, ThresholdG, 0xFF, ThresholdType.Binary); //matbuffer3= matbuffer2=matbuffer0>ThresholdG? 0xFF:0 
                        matbuffer2.ConvertTo(matbuffer3, DepthType.Cv8U);
                        CvInvoke.BitwiseAnd(localmaxima, matbuffer3, matbuffer4); //matbuffer4= localmaxima & matbuffer3
                        //win1 = "final spots";
                        //Program.show_grayimage(matbuffer3, win1, Nrows, Ncols);

                        flag_record = false;
                        NFid[nslice] = 0; //generate a fresh list of markers 
                        arr_matbuffer4 = matbuffer4.GetData();  //from mat to array, stores the fiducials according to general procedure
                        for (int row = marginy + exclude_radius; row < Nrows - marginy - exclude_radius; row++)
                        {
                            for (int col = marginx + exclude_radius; col < Ncols - marginx - exclude_radius; col++)
                            {
                                if (Convert.ToByte(arr_matbuffer4.GetValue(row, col)) == 0xFF && NFid[nslice] < NfidMax)
                                {
                                    flag_record = true;
                                    for (int ind_exc = NFid[nslice] - 1; ind_exc >= 0; ind_exc--)
                                    {
                                        if (Math.Abs(row - locations[nslice, ind_exc].row) <= exclude_radius && Math.Abs(col - locations[nslice, ind_exc].col) <= exclude_radius)
                                        {
                                            if (Math.Sqrt(Math.Pow((row - locations[nslice, ind_exc].row), 2) + Math.Pow((col - locations[nslice, ind_exc].col), 2)) <= exclude_radius)
                                            {
                                                flag_record = false;
                                                break;
                                            }
                                        }
                                    }
                                    if (flag_record)
                                    {
                                        locations[nslice, NFid[nslice]] = new tp(row, col);// write in table the location of fiducial found (y,x)
                                        NFid[nslice]= NFid[nslice]+1;
                                    }
                                }
                            }
                        }
                        //level_expected = Math.Min(1, (2 * cluster_size / Nrows) ^ 2) * NFid[ncenter];

                    }//if optical_test
                    win1 = "Fiducials Marked";
                    Program.show_circled_image(slice_mat, win1, Nrows, Ncols, locations, nslice, NFid[nslice], attention_size, fid_count, fidx, fidy,fidn, ref smatch, NfidMax);
                }
                CvInvoke.DestroyWindow(win1);

                Console.WriteLine("Processing cluster matching ...");

                //find approx. z height of fiducials according to center slice and neighbor slices to improve cluster tracking
                if (IterationNum == 0)// && ClusterAlign.Window1.loadfidFile == "")
                {
                    //z_values = nogaps.estimate_z(xisRotation, ncenter, Ncols, Nrows, locations, NFid, PreAlignmentTolx, PreAlignmentToly, ref tiltangles);
                    z_values = new double[NFid[ncenter]];
                }


                //In a second iteration, erase confusing points in table locations that do not fit to rigid body model according to previous iteration
                 double found_distance_sqr = Math.Pow(0.05 * Ncols, 2);


                // **** Determine structure: (1) collect vectors between nearby fiducials to be coupled ****
                float dist_points;
                for (nslice = 0; nslice < Nslices; nslice++)
                {
                    xc_ns = (int)(nslice == ncenter ? xc_ncenter : xc);
                    yc_ns = (int)(nslice == ncenter ? yc_ncenter : yc);
                    Nsvector[nslice] = 0;
                    for (int point1 = 0; point1 < NFid[nslice]; point1++)
                    {
                        col1 = (locations[nslice, point1].col - xc_ns) * factordx(nslice, tiltangles);  //factoring according to 1/cos(Theta) helps stabilize the vectors and helps decrease search window for matching
                        row1 = (locations[nslice, point1].row - yc_ns) * factordy(nslice, tiltangles);
                        for (int point2 = 0; point2 < NFid[nslice]; point2++)
                        {
                            col2 = (locations[nslice, point2].col - xc_ns) * factordx(nslice, tiltangles);
                            row2 = (locations[nslice, point2].row - yc_ns) * factordy(nslice, tiltangles);
                            dx_norm = (col2 - col1);
                            dy_norm = (row2 - row1);
                            dist_points = MathF.Sqrt(MathF.Pow(col1 - col2, 2) + MathF.Pow(row1 - row2, 2));
                            if (dist_points < (float)cluster_size && dist_points > 3 * fidsize && point1 != point2)
                            {
                                //Addition of estimated image shift Dx, Dy allows decreasing matching tolerance
                                svectors[nslice, Nsvector[nslice]] = new svector(point1, point2, (int)MathF.Round(dx_norm), (int)MathF.Round(dy_norm), (int)MathF.Round(col1 + (float)(Dx_vect[nslice] * factordx(nslice, tiltangles))), (int)MathF.Round(row1 + (float)(Dy_vect[nslice] * factordy(nslice, tiltangles))));
                                svectors_radius[nslice, Nsvector[nslice]] = (int)MathF.Round(0.5f*dist_points/tolerance_radius); //Sort key to accelerate search
                                Nsvector[nslice]++;
                            }
                        }
                    }
                }
                // **** Determine structure: (2) match between vectors across slices ****
                //smatch Initialize to -1 value
                for (int nt1 = 0; nt1 < Nslices; nt1++)
                {
                    for (int nt2 = 0; nt2 < NfidMax; nt2++)
                    {
                        smatch[nt1, nt2] = -1; ;
                    }
                }


                //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

                int nid;
            if (!isfidRequest)
            {
                    //Track fiducials, acquire results in table smatch
                    int Max_no_of_threads = 19;
                    int[] th_stop = new int[Max_no_of_threads];
                    Array.Clear(th_stop, 0, Max_no_of_threads);
                    Thread[] objThread = new Thread[Max_no_of_threads];
                    // Select fiducial in the center slice as center of cluster
                    for (int n1 = 0; n1 < NFid[ncenter]; n1++)
                    {
                        int th_idx = n1 % Max_no_of_threads;
                        objThread[th_idx] = new Thread(() => tracking(th_idx,n1, ncenter, NfidMax, NFid[ncenter], Nslices, NumberofLevels, (int)((N_minimum_tracked_fiducials_percent / 100.0) * Nslices), ref Nsvector, ref svectors, ref svectors_radius, ref smatch, ref th_stop));
                        objThread[th_idx].Start();

                        if (th_idx == Max_no_of_threads - 1 || n1 == NFid[ncenter]-1)
                        {
                            // Wait until thread is finished.
                            int count_finished_threads = 0;
                            while(count_finished_threads < th_idx+1)
                            {
                                count_finished_threads = 0;
                                for (int i = 0; i <= th_idx; i++)
                                {
                                    count_finished_threads= count_finished_threads+th_stop[i];
                                }

                            }
                            Array.Clear(th_stop, 0, Max_no_of_threads);
                            Console.Write("{0}%  ",(int)(100.0*n1/(NFid[ncenter] - 1)));
                            if (n1 == NFid[ncenter] - 1)
                            { Console.Write("\n"); }
                        }

                    }  //for n1

                    //Reduce too many fiducials with insufficient tracking
                    int count_collisions = 0;
                    int[] track_counter = new int[NFid[ncenter]];
                    Array.Clear(track_counter, 0, NFid[ncenter]);
                    for (int n1 = 0; n1 < NFid[ncenter]; n1++)
                    {
                        for (nslice = 0; nslice < Nslices; nslice++)  
                        {
                            count_collisions = 0;
                            for (nid = 0; nid < NFid[nslice]; nid++)
                            {
                                if (smatch[nslice, nid] == n1)
                                {
                                    track_counter[n1] = track_counter[n1] +1;
                                    count_collisions++; //only one in a slice is expected. Don't add weight to wrong matching.
                                }

                            }
                            if (count_collisions>1 && !ignore_collisions) //redundant, see tracking function where better collision blockage is implemented
                            {
                                track_counter[n1] = track_counter[n1] - count_collisions;
                                for (nid = 0; nid < NFid[nslice]; nid++)
                                {
                                    if (smatch[nslice, nid] == n1)
                                    {
                                        smatch[nslice, nid] = -1; //cancel all matches in this slice with this fiducial id to avoid hopping between markers
                                    }

                                }

                            }
                        }

                    }
                    int count100 = 0;
                    int[] countper = new int[Nslices];
                    Array.Clear(countper, 0, Nslices);
                    for (int n1 = 0; n1 < NFid[ncenter]; n1++)
                    {
                        if (track_counter[n1] == Nslices)
                        {
                            count100++;
                        }
                        else
                        {
                            countper[track_counter[n1]]++; //histogram of filling counts for all trackings
                        }    
                    }
                    if (IterationNum > 0) //In second interation, after imrpoved optical tracking 
                    {
                        if (count100 >= 5)
                        {
                            N_minimum_tracked_fiducials_percent = 100; //despite user request, if there are enough fiducials for 100% tracking modify threshold to 100%
                            Console.WriteLine("Now using 100% optical tracking threshold");
                        }
                        else
                        {
                            int countsum = 0;
                            for (int n2=Nslices-1; n2>= (int)((N_minimum_tracked_fiducials_percent / 100.0) * Nslices); n2--)
                            {
                                countsum += countper[n2];
                                if (countsum>=9) //Best tracking of 9 fiducials is assumed enough, so tracking threshold will increase accordingly
                                {
                                    N_minimum_tracked_fiducials_percent = (int)(100.0*n2/Nslices);
                                    Console.WriteLine("Now using {0}% optical tracking threshold", N_minimum_tracked_fiducials_percent);
                                    break;
                                }
                            }
                        }
                    }
                    for (int n1 = 0; n1 < NFid[ncenter]; n1++)
                    {
                        if (track_counter[n1] < (int)((N_minimum_tracked_fiducials_percent / 100.0) * Nslices))
                        {
                            for (nslice = 0; nslice < Nslices; nslice++)  //all other slices
                            {
                                for (nid = 0; nid < NFid[nslice]; nid++)
                                {
                                    if (smatch[nslice, nid] == n1)
                                    {
                                          smatch[nslice, nid] = -1; 
                                    }

                                }
                            }
                        }

                    }
                     //cancelled: exclude unproductive slices assuming they are on high tilt angles
                    minslice = 0;
                maxslice = Nslices - 1;
                //find preliminary range of slices
                for (nslice = 0; nslice < ncenter; nslice++)
                {
                    count = 0;
                    for (int n = 0; n < NfidMax; n++)
                    {
                        if (smatch[nslice, n] >= 0) { count++; }
                    }
                    if (count < 1 && nslice == minslice) { minslice++; }
                    else
                    { break; }
                }
                for (nslice = Nslices - 1; nslice > ncenter; nslice--)
                {
                    count = 0;
                    for (int n = 0; n < NfidMax; n++)
                    {
                        if (smatch[nslice, n] >= 0) { count++; }
                    }
                    if (count < 1 && nslice == maxslice) { maxslice--; }
                    else
                    { break; }
                }
                Console.WriteLine("Productive slices: " + minslice.ToString() + " to " + maxslice.ToString());
                minslice = 0; //overwrite minimum and maximum to get all possible points
                maxslice = Nslices - 1;

                // report locations of individual fiducials across all slices in fidx/y table
                //here match_counter is a table of number of matches for each fiducial id
                for (int n = 0; n < NfidMax; n++)
                {
                    match_counter[n] = 0;
                }
                for (int n = 0; n < NfidMax; n++)
                {
                    for (nslice = 0; nslice < Nslices; nslice++)
                    {
                        if (smatch[nslice, n] >= 0)
                        { match_counter[smatch[nslice, n]]++; } // smatch[nslice, n] holds the fiducial id (identical to the index n in smatch[ncenter,n] and to the values in 
                    }
                }
                fid_count = 0;
                missing_location = false;
                for (int n = 0; n < NFid[ncenter]; n++) //loop over fiducial ids
                {
                    if (smatch[ncenter, n] >= 0) //if enough matching was found for this fiducial in ncenter slice
                    {
                            for (nslice = 0; nslice < Nslices; nslice++)
                            {
                                missing_location = true;
                                //nslice_neighbor = nslice + (nslice > ncenter ? -1 : 1);
                                for (int nn = 0; nn < NfidMax; nn++) //loop over local particle index in nslice frame (the order which particles were found, not the order of fiducial id)
                                {
                                    if (smatch[nslice, nn] == smatch[ncenter, n])  //if matching for this fiducial was found and locked (the value of smatch is the fiducial id)
                                    {
                                        fidx[nslice, fid_count] = locations[nslice, nn].col;
                                        fidy[nslice, fid_count] = locations[nslice, nn].row;
                                        fidn[nslice, fid_count] = n; //fiducial id (fid number in ncenter slice)
                                        missing_location = false;
                                    }
                                }
                                if (missing_location)
                                {
                                    //complete missing locations according to known vectors dx,dy in ncenter slice
                                    fidx[nslice, fid_count] = 0;
                                    fidy[nslice, fid_count] = 0;
                                    fidn[nslice, fid_count] = -1;
                                    /*min_distance = (int)Math.Floor(tolfactor_stability * Ncols);
                                    temp_distance = min_distance + 1;
                                    for (int ntag1 = 0; ntag1 < Nsvector[ncenter]; ntag1++) //search over all vectors in center slice
                                    {
                                        if (svectors[ncenter, ntag1].arb_id2 == n) //n is the fiducial id of the missing point in nslice slice, looking for vector in center slice holding this id as target
                                        {
                                            bool missing_location_little = true;
                                            for (int nn = 0; nn < NfidMax && missing_location_little; nn++) //loop over local fiducial indices in nslice frame, looking for point matching with souce point in considered vector
                                            {
                                                if (smatch[nslice, nn] == svectors[ncenter, ntag1].arb_id1) //if the matching fiducial id is the source point of the vector then fill the gap
                                                {
                                                    //look only for points proven to participate in the same cluster with n (so the structure is sure to conserve)
                                                    for (int ii=1; ii<=friend[n,0]; ii++)
                                                    {
                                                        if (friend[n,ii]== svectors[ncenter, ntag1].arb_id1)
                                                        {
                                                            missing_location_little = false;
                                                        }
                                                    }
                                                    if (!missing_location_little)
                                                    {
                                                        tempx = locations[nslice, nn].col + (int)MathF.Round(svectors[ncenter, ntag1].arrow_dx / factordx(nslice, tiltangles));
                                                        tempy = locations[nslice, nn].row + (int)MathF.Round(svectors[ncenter, ntag1].arrow_dy / factordy(nslice, tiltangles));
                                                        fidx[nslice, fid_count] = tempx;
                                                        fidy[nslice, fid_count] = tempy;
                                                    }
                                                }
                                            }
                                        }
                                    }*/

                                }
                            }
                            fid_count++;
                    }
                }

            } //end of if !isfidRequest

            minslice = 0;
            maxslice = Nslices - 1;


                //produce *.fid.txt with all observed that can be read by tomoalign
                writeIMODfidModel(fidx, fidy, Ncols, Nrows, minslice, maxslice, fid_count, FidFileName);

                Bfinal = new double[Nslices* fid_count, 5];
                //fill all gaps in fidx, fidy according to rigid body model
                double model_error=-1;
                double min_model_error = -1;
                for (int n = 0; n < 1; n++)
                {
                    model_error = nogaps.fillgaps(0.0, xisRotation ? 0.5 * Math.PI : 0, tiltangles, Nslices, fid_count, ncenter, Ncols, Nrows, ref fidx, ref fidy, ref fidn, ref Bfinal, ref Dx_vect, ref Dy_vect, locations, NFid, IterationNum, ref CorrectAngle, PreAlignmentTolx, PreAlignmentToly, ref xc_ncenter, ref yc_ncenter, ref report_phi, ref report_psi, basefilename);
                    Console.WriteLine("(Fit err {0:0.00} pixels)", model_error);
                    if (min_model_error < 0)
                    {
                        min_model_error = model_error;
                    }
                    else if (model_error < min_model_error)
                    {
                        min_model_error = model_error;
                        break;
                    }
                }
                reportfobj.WriteLine("Iteration:{0}", IterationNum);
                reportfobj.WriteLine("phi={0:0.0}   psi={1:0.0}  [deg]", report_phi*180/Math.PI, report_psi * 180 / Math.PI);
                Console.WriteLine("Based on N={0} and optical detection in {1}% of the slices or higher:", NumberofLevels+1, N_minimum_tracked_fiducials_percent);
                reportfobj.WriteLine("Based on N={0} and optical detection in {1}% of the slices or higher:", NumberofLevels + 1, N_minimum_tracked_fiducials_percent);
                Console.WriteLine("Rigid body fitting error {0:0.0} pixles, tracking {1} fiducials.", model_error,fid_count);
                reportfobj.WriteLine("Rigid body fitting error {0:0.0} pixles, tracking {1} fiducials.", model_error, fid_count);
                if (fid_count>1)
                {
                    Console.WriteLine("Expected alignment residue for rigid body: {0:0.0} pixels", model_error / Math.Sqrt(fid_count));
                }
                writeNogapModel(ref Bfinal, fid_count, Nslices, NogapsFidFileName);//add table location to compare and fit with simulated points to extract Dx,Dy when tracked fiducials are missing

                if (cluster_visualize)
                {
                    for (nslice=0; nslice<Nslices; nslice++)
                    {
                        slices[nslice].ConvertTo(slice_mat, DepthType.Cv32F, 1, 0);
                        if (isMRCfile)
                        {
                            CvInvoke.Transpose(slice_mat, slice_mat); //switch x and y axis so afterwards x(col) is along 3dmod x and y[row] is along 3dmod y. The origin is upperleft corner here and is lowerleft corner in 3dmod but it is the same relative to the image
                        }
                        else
                        {
                            CvInvoke.Flip(slice_mat, slice_mat, FlipType.Vertical);
                        }

                        win1 = "Visualize one cluster";
                        Program.show_circled_image(slice_mat, win1, Nrows, Ncols, locations, nslice, NFid[nslice], attention_size, fid_count, fidx, fidy, fidn, ref smatch, NfidMax);
                        CvInvoke.WaitKey(5000);  //Wait for the key pressing event
                    }
                }

                attention_size = Math.Max((int)(model_error * Math.Sqrt(2) + fidsize),20);
                

                //Forcefill is the flag to use iterations. The knowledge acquired in iterations is the attention table.
                if (IterationNum==0 &&  NumofIterations > 1 && fid_count > 0 )
                {
                    Console.WriteLine("Building attention map and trying enhanced optical tracking.");
                    optical_test = true;
                    PseudoAttentionRequest = true;
                    //no, keep it! isAttentionRequest = false;
                    attention_NfidMax =fid_count;
                    //prepare Attention images based on best fiducials found so far, to be used in next iteration
                    for (nslice = 0; nslice < Nslices; nslice++)
                    {
                        PsAttentionslices[nslice] = new Mat(Nrows, Ncols, DepthType.Cv8U, 1);
                        PsAttentionslices[nslice].SetTo(new MCvScalar(255)); //fill all white
                        for (int n = 0; n < fid_count; n++)
                        {
                            lastpx = fidx[nslice, n];
                            lastpy = fidy[nslice, n];
                            CvInvoke.Circle(PsAttentionslices[nslice], new Point(lastpx, lastpy), attention_size, new MCvScalar(0), -1); //place filled circles with color 0
                        }
                    }
                }
                else
                { 
                    optical_test = false;
                    PseudoAttentionRequest = false;
                }



            } //end of for IterationNum
            string runmatlab = "";
            if (!Settings4ClusterAlign2.Default.export_normalali)
                runmatlab = "clusteralign_astra_reconstruct(" + (coswindow ? 1 : 0).ToString() + "," + Math.Round(report_phi * 180 / Math.PI).ToString() + "," + Math.Round(report_psi * 180 / Math.PI).ToString() + ",'" + out_filename + "','" + rawtltFilename +"','"+basefilename+ ".fit_err_by_slice.txt')";
            else
                runmatlab = "clusteralign_astra_reconstruct(" + (coswindow ? 1 : 0).ToString() + "," + 0.ToString() + "," + Math.Round(report_psi * 180 / Math.PI).ToString() + ",'" + normout_filename + "','" + rawtltFilename + "','" + basefilename + ".fit_err_by_slice.txt')";
            reportfobj.WriteLine("Command for reconstruction in Matlab:");
            reportfobj.WriteLine(runmatlab);

            reportfobj.Close();

            if (fid_count == 0)
            {
                Console.WriteLine("### Attempt failed. Suggestion: Try increasing the max number of fiducials or decreasing the tracking threshold. ###");
                //System.Environment.Exit(0);
            }
            else
            {
                Console.WriteLine("Range of slices: " + minslice.ToString() + " to " + maxslice.ToString());
                //generate ali file
                MrcStack myMrcStack = new MrcStack();
                if (!Settings4ClusterAlign2.Default.export_normalali)
                {
                    Console.WriteLine("Saving ali file (_jali.mrc) ..");
                    myMrcStack.Fexport(FileName, out_filename, ref Dx_vect, ref Dy_vect, ref report_phi);
                }
                else
                {
                    Console.WriteLine("Saving normal ali file for 3rd party software ..");
                    myMrcStack.Fexport(FileName, normout_filename, ref Dx_vect, ref Dy_vect, ref report_phi);
                }

                //generate reconstruction if requested by user
                if (ClusterAlign.Settings4ClusterAlign2.Default.add_reconst)
                {
                    try
                    {
                        #if Windows
                        Console.WriteLine("Reconstruction ..");
                        MLApp.MLApp matlab = new MLApp.MLApp();
                        //Console.WriteLine(@runmatlab);
                        matlab.Execute(@runmatlab);
                        #else
                        Console.WriteLine("For reconstruction past the command in Matlab: "+runmatlab);
                        #endif
                    }
                    catch
                    {
                        Console.WriteLine("Connection to Matlab failed, use Matlab directly (see output.txt)");
                    }
                }
             }

            Console.WriteLine("Program ended");

        }
 
    
    
    
        public class tp
        {
            public int row = 0;
            public int col = 0;
            public tp(int r, int c)
            {
                row = r;
                col = c;
            }
        }
        public class svector
        {
            public int arb_id1;
            public int arb_id2;
            public int arrow_dx;
            public int arrow_dy;
            public int origin_x;
            public int origin_y;
            public svector(int a,int b, int c, int d, int e, int f)
            {
                arb_id1 = a; arb_id2 = b; arrow_dx = c; arrow_dy = d; origin_x = e; origin_y = f;
            }
        }
        public class permutation
        {
            public int p1;
            public int p2;
            public int p3;
            public int p4;
            public permutation(int a,int b,int c,int d)
            {
                p1 = a; p2 = b; p3 = c; p4 = d;
            }
        }
        
        public static bool IsSvMatch(svector sv1, svector sv2, double delta_z, double theta1, double theta2,int ns) //Older version, not used! see IsSvMatch_faster below
        {
            //one of the vectors should be in ncenter slice, so delta_z should be known and tentavely assumed the same for the other vector in another slice
            //arrows dx,dy are already compensated for contraction by cos(theta), overall for y axis roation x'=cos(theta)*x+sin(theta)*z should conserve between slices
            //used with: IsSvMatch(svec - ref in ncenter slice, svectors[nslice, ntag], z_values[svec.arb_id2]- z_values[svec.arb_id1], tiltangles[ncenter], tiltangles[nslice]))
            //deltaz is in the same order of the arrowx/y: point2-point1.
            if (Math.Abs(sv1.origin_x - sv2.origin_x) > match_tolerance[ns].origin_x || Math.Abs(sv1.origin_y - sv2.origin_y) > match_tolerance[ns].origin_y)
            {
                return false;
            }
            int rot_delx = 0;
            int rot_dely = 0;
            if (xisRotation) { rot_dely = (int)(delta_z * (Math.Sin(theta1) - Math.Sin(theta2))); }
            else { rot_delx = (int)(delta_z * (Math.Sin(theta1) - Math.Sin(theta2))); }
            if (Math.Abs(sv1.arrow_dx - sv2.arrow_dx + rot_delx) < match_tolerance[ns].arrow_dx && Math.Abs(sv1.arrow_dy - sv2.arrow_dy + rot_dely) < match_tolerance[ns].arrow_dy)
            {
                //if (Math.Abs(sv1.origin_x - sv2.origin_x) < match_tolerance.origin_x && Math.Abs(sv1.origin_y - sv2.origin_y) < match_tolerance.origin_y)
                //{
                //    return true;
                //}
                return true; 
                //else
                //{
                //    return false;
                //}
            }
            else
            {   return false;
            }
        }
        public static bool IsSvMatch_faster(svector sv1, svector sv2,int ns)
        {
            if (Math.Abs(sv1.origin_x - sv2.origin_x) > match_tolerance[ns].origin_x || Math.Abs(sv1.origin_y - sv2.origin_y) > match_tolerance[ns].origin_y)
            {
                return false;
            }
            if (Math.Abs(sv1.arrow_dx - sv2.arrow_dx ) >= match_tolerance[ns].arrow_dx || Math.Abs(sv1.arrow_dy - sv2.arrow_dy) >= match_tolerance[ns].arrow_dy)
            {
                return false;
            }
            else
            {
                return true;
            }
        }
        public static bool IsSvNear(svector sv1, svector sv2)
        {
            if (Math.Abs(sv1.origin_x - sv2.origin_x) < cluster_size && Math.Abs(sv1.origin_y - sv2.origin_y) < cluster_size)
            { return true; }
            else
            { return false; }
        }
        public static float Distance(int nslice1, int n1, int nslice2, int n2, ref tp[,] locations)
        {
            return MathF.Sqrt(MathF.Pow(locations[nslice1, n1].col - locations[nslice2, n2].col, 2) + MathF.Pow(locations[nslice1, n1].row - locations[nslice2, n2].row, 2));
        }
        public static float factordx(int nslice, double[] tiltangles)
        {
            if (coswindow)
            { return 1F; }
            else
            {
                if (xisRotation)
                { return 1F; }
                else
                { return (1F / MathF.Cos((float)tiltangles[nslice])); }
            }
        }
        public static float factordy(int nslice, double[] tiltangles)
        {
            if (coswindow)
            { return 1F; }
            else
            {
                if (xisRotation)
                { return (1F / MathF.Cos((float)tiltangles[nslice])); }
                else
                { return 1F; }
            }
        }
        static System.Drawing.Point mypoint = new System.Drawing.Point((int)0, (int)0);


        public static void tracking(int i, int n1, int ncenter, int NfidMax, int Nfid_center, int Nslices, int NumberofLevels, int min_slice_count, ref int[] Nsvector, ref svector[,] svectors, ref int[,] svectors_radius, ref int[,] smatch, ref int[] th_stop)
        {
            //The purpose is to fill up smatch table in many threads
            int tempNumberofLevels, startNumberofLevels, endNumberofLevels;
            if (NumberofLevels<=0) //auto =-1
            {
                startNumberofLevels = 5;
                endNumberofLevels = (int)(Nfid_center-1);
            }
            else //determined by user
            {
                startNumberofLevels = NumberofLevels;
                endNumberofLevels = NumberofLevels;
            }

            int list_size=0;
            int[] cluster_candidate_list = new int[Nfid_center];
            int[,] match_score = new int[Nslices,NfidMax];
            //Array.Clear(match_score, 0, NfidMax*Nslices);
            int nid;
            for (int vector_num = 0; vector_num < Nsvector[ncenter]; vector_num++)
            {
                if (svectors[ncenter, vector_num].arb_id1 == n1)
                {
                    cluster_candidate_list[list_size] = vector_num;
                    list_size++;
                }
            }
            int track_counter = 1;
            int s_radius = 0;
            for (int nslice = 0; nslice < Nslices; nslice++)  //all other slices
            {
                if (nslice != ncenter)
                {
                    int vis_ind = 0;
                    for (int list_num1 = 0; list_num1 < list_size; list_num1++)
                    {
                        svector svec = svectors[ncenter, cluster_candidate_list[list_num1]];
                        s_radius = svectors_radius[ncenter, cluster_candidate_list[list_num1]];
                        for (int vector_num2 = 0; vector_num2 < Nsvector[nslice]; vector_num2++)
                        {
                            //if(IsSvMatch(svec, svectors[nslice, vector_num2], z_values[svec.arb_id2] - z_values[svec.arb_id1], tiltangles[ncenter], tiltangles[nslice]))
                            if (svectors_radius[nslice, vector_num2] == s_radius) //accelerate filtering unmatched cases
                            {
                                if (IsSvMatch_faster(svec, svectors[nslice, vector_num2], nslice))
                                {
                                    nid = svectors[nslice, vector_num2].arb_id1;
                                    match_score[nslice, nid] = match_score[nslice, nid] + 1;
                                    //VISUALIZE
                                    if (cluster_visualize && n1 == n1visualize) 
                                    {
                                        visualize[nslice, vis_ind, 0] = svectors[nslice, vector_num2].arb_id1;//nid
                                        visualize[nslice, vis_ind, 1] = svectors[nslice, vector_num2].arb_id2;
                                        vis_ind++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            int best_n2 = -1;
            int score = 0;
            bool check_collision = false;
            bool collision = false;
            int best_tempNumberofLevels = 5;
            int score_tempNumberofLevels = 0;
            for (tempNumberofLevels = startNumberofLevels; tempNumberofLevels <= endNumberofLevels; tempNumberofLevels++)
            {
                track_counter = 0;
                for (int nslice = 0; nslice < Nslices; nslice++)  //all other slices
                {
                    if (nslice != ncenter)
                    {
                        score = 0;
                        check_collision = false;
                        collision = false;
                        for (nid = 0; nid < NfidMax; nid++)
                        {
                            if (match_score[nslice, nid] > tempNumberofLevels)
                            {
                                if (!check_collision)
                                {
                                    score = match_score[nslice, nid];
                                    if (!ignore_collisions) check_collision = true;
                                }
                                else
                                {
                                    collision = true;
                                }
                            }
                        }
                        if (score > tempNumberofLevels && !collision)
                        {
                            track_counter++;
                        }
                    } // if (nslice != ncenter)
                }
                if (track_counter>= score_tempNumberofLevels && track_counter>=min_slice_count) //try to increase NumberofLevels (confidence in tracking) but achieve minimum tracking count 
                {
                    best_tempNumberofLevels = tempNumberofLevels;
                    score_tempNumberofLevels = track_counter;
                }
            }
            if (NumberofLevels <= 0) //auto =-1
            {
                tempNumberofLevels = best_tempNumberofLevels;
            }
            else
            {
                tempNumberofLevels = NumberofLevels;
            }
            for (int nslice = 0; nslice < Nslices; nslice++)  //all other slices
            {
                if (nslice != ncenter)
                {
                    score = 0;
                    best_n2 = -1;
                    check_collision = false;
                    collision = false;
                    for (nid = 0; nid < NfidMax; nid++)
                    {
                        if (match_score[nslice, nid] > tempNumberofLevels)
                        {
                            if (!check_collision)
                            {
                                score = match_score[nslice, nid];
                                best_n2 = nid;
                                if (!ignore_collisions) check_collision = true;
                            }
                            else
                            {
                                collision = true;
                            }
                        }
                    }
                    if (score > tempNumberofLevels && !collision)
                    {
                        lock (smatch)
                        { smatch[nslice, best_n2] = n1; }
                        // smatch_strength1[nslice, best_n2] = best_score;
                        track_counter++;
                    }
                } // if (nslice != ncenter)
                else
                {
                    lock (smatch)
                    { smatch[ncenter, n1] = n1; }
                }
            }
            lock (th_stop)
            { th_stop[i]= 1; }

        }

        //Older, recursive code to match clusters, was less effective.
        /*public static bool SeekLevel(int level, ref int first_ntag, int last_ntag, int[] ntag_array, int Nsvectors_nslice, int ncenter, ref int NumberofLevels, ref int n1, ref int nslice, ref int inmaster, ref svector[,] svec_perm)
        {
            bool passed;
            svector svec = svec_perm[inmaster, level - 1];
            if (svec == null) return false;
            for (int ntag = 0; ntag < Nsvectors_nslice; ntag++)
            {
                bool wasbefore = false;
                for (int j = 1; j < level; j++) //compare to previous ntags
                {
                    if (ntag == ntag_array[j]) wasbefore = true;
                }
                if (!wasbefore && (level==1 || ( svectors[nslice, last_ntag].arb_id2 == svectors[nslice, ntag].arb_id1 && IsSvNear(svectors[nslice, first_ntag], svectors[nslice, ntag] ))))
                {
                    if (IsSvMatch(svec, svectors[nslice, ntag], z_values[svec.arb_id2]- z_values[svec.arb_id1], tiltangles[ncenter], tiltangles[nslice],nslice)) 
                    {
                        if (level == 1) first_ntag = ntag;
                        last_ntag = ntag;
                        ntag_array[level] = ntag;
                        if (level < NumberofLevels)
                        {
                            passed=SeekLevel(level + 1, ref first_ntag, last_ntag, ntag_array, Nsvectors_nslice, ncenter, ref NumberofLevels,ref n1, ref nslice, ref inmaster, ref svec_perm);
                            if (passed) return true;//tell the previous levels to stop searcing since a match is found
                        }
                        else
                        {
                             //nextsv[nslice] = svectors[nslice, ntag].arb_id2;
                            return true;
                        }
                    }
                }
            }
            return false;
        }*/



        public static void show_grayimage(Mat source, string win1, int Nrows, int Ncols)
        {
            Mat showfield = new Mat(Nrows, Ncols, DepthType.Cv8U, 1);
            double minVal = 0;
            double maxVal = 0;
            System.Drawing.Point locationh = new System.Drawing.Point(0, 0);
            System.Drawing.Point locationl = new System.Drawing.Point(0, 0);
            CvInvoke.MinMaxLoc(source, ref minVal, ref maxVal, ref locationl, ref locationh);
            double scale = 255 / (maxVal - minVal);
            CvInvoke.ConvertScaleAbs(source, showfield, scale, -scale * minVal);
            //CvInvoke.EqualizeHist(source, showfield);
            CvInvoke.NamedWindow(win1, WindowFlags.KeepRatio); //Create the window using the specific name
            CvInvoke.Imshow(win1, showfield); //Show the image
            CvInvoke.WaitKey(500);  //Wait for the key pressing event
            CvInvoke.DestroyWindow(win1);
        }
        public static void show_circled_image(Mat source, string win1, int Nrows, int Ncols, tp[,] locations, int nslice, int Nfid, int attention_size, int fid_count, int[,] fidx, int[,] fidy, int[,] fidn, ref int[,] smatch, int NfidMax)
        {
            Mat showfieldV0 = new Mat(Nrows, Ncols, DepthType.Cv8U, 1);  //1 channel to be scaled
            Mat showfield = new Mat(Nrows, Ncols, DepthType.Cv8U, 3);  //3 channels (RGB)
            Mat shsource = new Mat(Nrows, Ncols, DepthType.Cv32F, 1);
            double minVal = 0;
            double maxVal = 0;
            findminmaxVal(ref source, ref minVal, ref maxVal, (int)((Nrows / 250)* (Nrows / 250)));
            double scale = 255 / (maxVal - minVal);
            CvInvoke.ConvertScaleAbs(source, showfieldV0, scale, -scale * minVal);
            CvInvoke.CvtColor(showfieldV0, showfield, ColorConversion.Gray2Bgr); //copy from one channel (gray) to 3 channels (color)
            MCvScalar myred = new MCvScalar(0, 0, 0xFF);
            MCvScalar mygreen = new MCvScalar(0, 0xFF, 0);
            MCvScalar mywhite = new MCvScalar(0xFF, 0xFF, 0xFF);
            for (int n=0; n< Nfid; n++) 
            {
                CvInvoke.Circle(showfield, new System.Drawing.Point(locations[nslice, n].col, locations[nslice, n].row), (int)(Nrows/250), myred, 1 + (int)(Nrows / 2500), LineType.AntiAlias);
            }
            int lastpx,lastpy;

            if (attention_size>0)
            {
                for (int n = 0; n < fid_count; n++)
                {
                    lastpx = fidx[nslice, n];
                    lastpy = fidy[nslice, n];
                    if (fidn[nslice, n]>=0 )
                    { CvInvoke.Circle(showfield, new System.Drawing.Point(lastpx, lastpy), attention_size, mygreen, 1+(int)(Nrows / 1500), LineType.AntiAlias); } //tracked
                    if (fidn[nslice, n]==-2 )
                    { CvInvoke.Circle(showfield, new System.Drawing.Point(lastpx, lastpy), attention_size, mywhite, 1 + (int)(Nrows / 1500), LineType.AntiAlias); }//extrapolated
                }
            }
            if (cluster_visualize)
            {
                for (int nd=0; nd<50000; nd++)
                {
                    if (smatch[nslice, visualize[nslice, nd, 0]] == n1visualize)
                    {
                        int p1x = locations[nslice, visualize[nslice, nd, 0]].col;
                        int p1y = locations[nslice, visualize[nslice, nd, 0]].row;
                        int p2x = locations[nslice, visualize[nslice, nd, 1]].col;
                        int p2y = locations[nslice, visualize[nslice, nd, 1]].row;
                        CvInvoke.Line(showfield, new Point(p1x, p1y), new Point(p2x, p2y), mygreen,2 , LineType.AntiAlias);
                    }
                }

            }

            CvInvoke.PutText(showfield,"Slice no="+nslice.ToString(),new Point((int)(Nrows/100),(int)(Ncols/50)),FontFace.HersheyPlain,(int)(Nrows/500),new MCvScalar(255,255,0),2+ (int)(Nrows / 1500));
            CvInvoke.NamedWindow(win1, WindowFlags.KeepRatio); // WindowFlags.Fullscreen   Create the window using the specific name
            CvInvoke.Imshow(win1, showfield); //Show the image
            CvInvoke.WaitKey(100);  //Wait for the key pressing event
            //CvInvoke.DestroyWindow(win1);
        }
        public static void SpatialGradient(Mat source, ref Mat gradx, ref Mat grady)
        {
            //showfield: prepare 8bits map stretching of the source, since SpatialGradient calcualtion only works with cv8U
            Mat showfield = new Mat(source.Rows, source.Cols, DepthType.Cv8U, 1);
            Mat temp_gradx = new Mat(source.Rows, source.Cols, DepthType.Cv16S, 1);
            Mat temp_grady = new Mat(source.Rows, source.Cols, DepthType.Cv16S, 1);
            double minVal = 0;
            double maxVal = 0;
            System.Drawing.Point locationh = new System.Drawing.Point(0, 0);
            System.Drawing.Point locationl = new System.Drawing.Point(0, 0);
            CvInvoke.MinMaxLoc(source, ref minVal, ref maxVal, ref locationl, ref locationh);
            double scale = 255 / (maxVal - minVal);
            CvInvoke.ConvertScaleAbs(source, showfield, scale, -scale * minVal);
            CvInvoke.SpatialGradient(showfield, temp_gradx, temp_grady, 3); //ksize=3
            temp_gradx.ConvertTo(gradx, DepthType.Cv32F);
            temp_grady.ConvertTo(grady, DepthType.Cv32F);
        }
        public static float FindThreshold(Mat source, IInputArray mask, int NfidMax, int kmask_area, bool bright_features)
        {
            double minVal = 0;
            double maxVal = 0;
            double mutiply_factor = 1;// mask==null? 1.2:1.2;// if mask is still not ready then prepare more candidates of fiducials, so expand their allowed numbers
            double approx_count;
            MCvScalar avgmat;
            Mat shsource = new Mat(source.Rows, source.Cols, DepthType.Cv32F, 1);
            System.Drawing.Point locationh = new System.Drawing.Point(0, 0);
            System.Drawing.Point locationl = new System.Drawing.Point(0, 0);
            findminmaxVal(ref source, ref minVal, ref maxVal, kmask_area);
            avgmat = CvInvoke.Mean(source);
            minVal = avgmat.V0;
            if (maxVal - minVal == 0)
            {
                Console.WriteLine("empty image");
                return 0; }
            double threshold_t = 0;
            Mat thresh = new Mat(source.Rows, source.Cols, DepthType.Cv8U, 1);
            source.ConvertTo(shsource, DepthType.Cv32F,255d/(maxVal-minVal),-minVal* 255d / (maxVal - minVal));
            shsource.ConvertTo(thresh, DepthType.Cv8U);
            threshold_t = CvInvoke.Threshold(thresh, thresh, 0, 255, ThresholdType.Triangle); //Finds threshold by Otsu/triangle algorithm
            avgmat = CvInvoke.Sum(thresh);
            approx_count = (avgmat.V0 / 255.0) / kmask_area;
            int counter = (approx_count > NfidMax && threshold_t<5) ?1:0;
            float threshold = 0;
            do
            {
                if (counter>0 )
                {
                    threshold_t = threshold_t+5;
                    shsource.ConvertTo(thresh, DepthType.Cv8U);
                    CvInvoke.Threshold(thresh, thresh, threshold_t, 255, ThresholdType.Binary);
                    if (mask != null) 
                    {
                        CvInvoke.Multiply(thresh, mask, thresh, 1.0, DepthType.Cv8U); 
                    }
                }
                threshold = (float)(threshold_t * (maxVal - minVal) / 255d + minVal);
                //display for debugging
                //string win2 = "check triangle threshold";
                //CvInvoke.NamedWindow(win2, WindowFlags.KeepRatio);
                //Program.show_grayimage(thresh, win2, source.Rows, source.Cols);
                //using SimpleBlobDetector class
                VectorOfVectorOfPoint contours = new VectorOfVectorOfPoint();
                Mat hierarchy = new Mat();
                CvInvoke.FindContours(thresh, contours, hierarchy, RetrType.External, ChainApproxMethod.ChainApproxSimple);
                counter = contours.Size; //note: count will be reasonable only if the white objects are enough sparse, otherwise could be low but the approx_count will be high
                //Console.WriteLine("Counted blobs=" + counter.ToString() + "  threshold=" + threshold.ToString());/////
            } while (counter > NfidMax);

            avgmat=CvInvoke.Sum(thresh);
            approx_count = (avgmat.V0 / 255.0) / kmask_area;
            if (counter>0 && approx_count< NfidMax*20) return threshold; //normally will exit here

            //As back plan count according to histogram
            int size = 256;
            int[] channels = new int[] { 1 };
            int[] histsize = new int[] { size };
            Mat hist = new Mat(1, size, DepthType.Cv16S, 1);
            Array arr_hist = new int[size,1];
            bool accumulate = false;
            //IInputArray mask = null; //null means to ignore
            float[] histrange=new float[] { (float)minVal, (float)maxVal };
            VectorOfMat sources=new VectorOfMat(source);
            sources.Push(source);

            CvInvoke.CalcHist(sources, channels, mask, hist, histsize, histrange , accumulate);
            arr_hist = hist.GetData();  //from mat to array
            int sum = 0;
             if (bright_features) {
                for (int ind = size - 1; ind >= 0; ind--)
                {
                    sum += Convert.ToInt32(arr_hist.GetValue(ind, 0))>= kmask_area ? Convert.ToInt32(arr_hist.GetValue(ind, 0)):0;
                    if (sum > (int)(NfidMax * kmask_area* mutiply_factor))// || (chunk<0.2*chunk_max && mask!=null))
                    {
                         break;
                    } else
                    {
                        threshold = (float)ind * (float)(maxVal - minVal) / (float)size + (float)minVal; //update before exceeds number
                    }
                }
            } 
            else
            {
                for (int ind = 0; ind <size; ind++)
                {
                    sum += Convert.ToInt32(arr_hist.GetValue(ind, 0))>= kmask_area ? Convert.ToInt32(arr_hist.GetValue(ind, 0)):0;
                    if (sum > (int)(NfidMax * kmask_area* mutiply_factor)) // || (chunk < 0.2 * chunk_max && mask!=null))
                    {
                        break;
                    }
                    else
                    {
                        threshold = (float)ind * (float)(maxVal - minVal) / (float)size + (float)minVal; //update before exceeds number
                    }
                }
            }
            Console.WriteLine("Alternative plan: threshold=" + threshold.ToString()); 
            return threshold;
        }

        public static void findminmaxVal(ref Mat source, ref double minVal, ref double  maxVal,int spot_area)
        {
            //Steps to find suitable mival/maxval for effectuive reduction to 8bit image
            int Nrow = source.Rows;
            int Ncol = source.Cols;
            int count=0;
            double minVal0=0, maxVal0=0; 
            System.Drawing.Point locationh = new System.Drawing.Point(0, 0);
            System.Drawing.Point locationl = new System.Drawing.Point(0, 0);
            Mat shsource = new Mat(Nrow, Ncol, DepthType.Cv32F, 1);
            //double dummy = 0;
            //Mat kern = CvInvoke.GetStructuringElement(ElementShape.Ellipse, new Size(3, 3), new Point(2, 2));
            //CvInvoke.Dilate(source, shsource, kern, anchor: new System.Drawing.Point(-1, -1), 1, BorderType.Replicate, new MCvScalar(0));
            //CvInvoke.MinMaxLoc(shsource, ref minVal0, ref dummy, ref locationl, ref locationh);
            //CvInvoke.Erode(source, shsource, kern, anchor: new System.Drawing.Point(-1, -1), 1, BorderType.Replicate, new MCvScalar(0));
            //CvInvoke.MinMaxLoc(shsource, ref dummy, ref maxVal0, ref locationl, ref locationh);
            CvInvoke.MinMaxLoc(source, ref minVal0, ref maxVal0, ref locationl, ref locationh);
            minVal = minVal0;
            maxVal = maxVal0;
            while (count <= Nrow)
            {
                minVal=minVal + (maxVal0 - minVal0) / 64;
                CvInvoke.Threshold(source, shsource, minVal, 0xFF, ThresholdType.BinaryInv); //shsource= source<val? 0xFF:0 
                count = CvInvoke.CountNonZero(shsource);
            }
            minVal = minVal - (maxVal0 - minVal0) / 64;
            count = 0;
            while (count <= spot_area*3)
            {
                maxVal = maxVal - (maxVal0 - minVal0) / 64;
                CvInvoke.Threshold(source, shsource, maxVal, 0xFF, ThresholdType.Binary); //shsource= source>val? 0xFF:0 
                count = CvInvoke.CountNonZero(shsource);
            } 
            maxVal = maxVal + (maxVal0 - minVal0) / 64;
        }

        public static void writeIMODfidModel(int[,] fidx, int[,] fidy, int width, int height, int minslice, int maxslice, int fid_count,string fidFileName)
        {
            System.IO.FileStream outfile = File.Create(fidFileName);
            var sr = new StreamWriter(outfile,System.Text.UTF8Encoding.Default);
            sr.NewLine = "\n";
            int count;
            for (int i = 0; i < fid_count; i++)
            {
                count = 0;
                for (int j = minslice; j <= maxslice; j++)
                {
                    if (fidx[j, i]!=0 || fidy[j, i]!=0)
                    { count++; }
                }
                //sr.WriteLine("contour " + i.ToString() + " 0 " + count.ToString());
                if (true)
                {
                    for (int j = minslice; j <= maxslice; j++)
                    {
                        if (fidx[j, i] >0 || fidy[j, i] > 0 )
                        {
                            sr.WriteLine(0.ToString() + "\t" + i.ToString() + "\t" + (fidx[j, i] + 1).ToString() + "\t" + (fidy[j, i] + 1).ToString() + "\t" + j.ToString());
                            Console.WriteLine(0.ToString() + "\t" + i.ToString() + "\t" + (fidx[j, i] + 1).ToString() + "\t" + (fidy[j, i] + 1).ToString() + "\t" + j.ToString());
                        }
                    }
                }
            }
            sr.Close();
            Console.WriteLine("Your file is ready: " + fidFileName);
        }

        public static void readIMODfidModel(ref int[,] fidx, ref int[,] fidy, ref int[,] fidn, ref int fid_count, string fidFileName)
        {
            int count = -1;
            StreamReader streamreader = new StreamReader(@fidFileName);
            char[] delimiter = new char[] { ' ' };// { '\t' };
            String temproot;
            string[] temp;
            int a, b, c, d,e;
            while (streamreader.Peek() > 0)
            {
                temproot = streamreader.ReadLine();
                temproot=System.Text.RegularExpressions.Regex.Replace(temproot, "\\s+", " ");
                temp =temproot.Split(delimiter);
                if (temp[0]=="")
                {
                    temp[0] = temp[1];
                    temp[1] = temp[2];
                    temp[2] = temp[3];
                    temp[3] = temp[4];
                    temp[4] = temp[5];
                }
                a = int.Parse(temp[0]);//object
                b = int.Parse(temp[1]);//contour / fiducial number (from 0)
                c = (int)float.Parse(temp[2]);//x
                d = (int)float.Parse(temp[3]);//y
                e = (int)float.Parse(temp[4]);//slice (from 0)
                if (b >= count) count=b+1;
                fidx[e, b] = c;
                fidy[e, b] = d;
                fidn[e, b] = b;
            }
            fid_count = count;
            streamreader.Close();
            Console.WriteLine("Continuing the fiducial table in file " + fidFileName);
            return;
         }

        public static void writeNogapModel(ref double[,] Bfinal,int fid_count, int Nslices, string NogapsfidFileName)
        {
            System.IO.FileStream outfile = File.Create(NogapsfidFileName);
            var sr = new StreamWriter(outfile, System.Text.UTF8Encoding.Default);
            sr.NewLine = "\n";
            for (int i = 0; i < fid_count*Nslices; i++)
            {
                sr.WriteLine(Bfinal[i,0].ToString() + "\t" + Bfinal[i, 1].ToString() + "\t" + Bfinal[i, 2].ToString() + "\t" + Bfinal[i, 3].ToString() + "\t" + Bfinal[i, 4].ToString());
            }
            sr.Close();
            Console.WriteLine("Your file is ready: " + NogapsfidFileName);
        }

    }
}
 